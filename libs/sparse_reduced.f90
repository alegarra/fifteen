!s Collection of predominantly sparse matrix structures and operations

module sparsem
use kinds
use hash
implicit none

integer::sparsem_msg=2	! message level. 0=no messages, 3=max messages
			! not yet implemented

!=====================================================================
! Type definitions
!=====================================================================

type densem     !traditional dense square matrix
   integer(i8) :: n
   integer(i8) :: filled    ! for compatibility mainly
   real (r8) ,pointer::x(:,:)=>null()
end type densem


type dense_symm       !upper stored symmetric dense matrix
   integer(i8) ::n
   integer(i8) :: filled    ! for compatibility mainly
   real (r8) ,pointer::x(:)=>null()
end type dense_symm


type sparse_hashm
    integer(i8):: n,&       ! for compatibility mainly
            nel,&       ! number of elements
            filled,&    ! number of filled elements
            status      ! 1 if ready to hash, 2 if in sorted order
     character:: storage_type ! storage type u upper(default), f full 
     real (rh) , pointer :: xreal(:,:)=>null() !real
     integer(i8),pointer :: xind(:,:)=>null() !int
end type sparse_hashm


type sparse_ija
     integer(i8) :: n,nel      !n is number of equations; nel is number of nonzeroes
     integer(i8), pointer::ia(:)=>null()      !will be ia(n+1)
     integer(i8), pointer::ja(:)=>null()      !will be ja(m)
     real (r8), pointer::a(:)=>null()     !will be a(m)
end type



type ll_element  ! basic element of linked-list
     integer(i8)::col
     real::val
     type (ll_element),pointer::next=>null()
end type    

type sparse_ll	! actual linked-list structure
    integer(i8) n
    type (ll_element),pointer::ia(:)=>null()
end type

!------------------------------------------
! Defaults and constants for subroutines
!------------------------------------------
! integer(i8)        ::     default_hash_size=5000,&
 integer(i8)        ::     default_hash_size=8192,& !AL
                       default_rounds = 1000  !maximum number of rounds
 real           ::     default_conv=1e-10,&   ! stopping criterion
                       default_relax=1.0,&      !relaxation factor
                       maxfill=.80,&          ! max fill of hash table
                       !hash_incr=1.5          ! increase of hash table
                       hash_incr=2          ! increase of hash table, AL
 logical        ::     default_zerosol=.true. ! zero solutions before 
                                              ! iterative solving	       
					      
					      					      
 integer(i8),parameter ::  conv_straight=0,&      ! straight conversion from hash
                                              !  to ija formats
                      &conv_upper_full=1      !	conversion from upper to full
		                              ! storage				      
 


!=====================================================================
!Interfaces
!=====================================================================

  interface init
      module procedure init_densem,&
                       init_dense_symm,&
                       init_sparse_hashm,&
                       init_sparse_ija
  end interface

  interface zerom
      module procedure zero_densem,&
                       zero_dense_symm,&
                       zero_sparse_hashm,&
                       zero_sparse_ija
  end interface

  interface reset
      module procedure deall_densem,&
                       deall_dense_symm,&
                       deall_sparse_hashm,&
                       deall_sparse_ija
  end interface

  interface addm
     module procedure add_densem, &
                      add_dense_symm,&
                      add_sparse_hashm
  end interface

  interface setm
     module procedure set_densem,&
                      set_dense_symm,&
                      set_sparse_hashm
  end interface

  interface getm
     module procedure get_densem,&
                      get_dense_symm,&
                      get_sparse_hashm,&
                      get_sparse_ija
  end interface


  interface	!external function
      function hashvr_old(dat,r,ind,m,n,mode,nr) !AL
      use kinds
      integer(i8) ::hashvr,r,mode,nr,m,n
      real (rh)::ind(:,:),dat(r)
      end function
  end interface


  contains

!=====================================================================
!Module subroutines
!=====================================================================

!----------------------------------------------------------------
!  initialize; 
!  necessary on systems where pointers not initialized at startup
!----------------------------------------------------------------
subroutine init_densem(x)
 ! Initializes a dense matrix
 type (densem) :: x
!
 x%n=0
 nullify(x%x)
end subroutine


subroutine init_dense_symm(x)
 ! Initializes dense symmetric matrix
 type (dense_symm) :: x
!
 x%n=0
 nullify(x%x)
end subroutine

subroutine init_sparse_hashm(x)
 ! Initializes sparse matrix in hash form
 type (sparse_hashm) :: x
!
 x%n=0; x%filled=0; x%status=0;x%storage_type='u'
 nullify(x%xind)
 nullify(x%xreal)
end subroutine

subroutine init_sparse_ija(x)
 ! Initializes sparse matrix in ija form
 type (sparse_ija) :: x
!
 x%n=0; x%nel=0
 nullify(x%ia,x%ja,x%a)
end subroutine

subroutine init_ll(x)
 ! Initializes sparse matrix in linked-list form
 type (sparse_ll) :: x
!
 x%n=0; nullify(x%ia)
end subroutine

!--------------
!  zero
!--------------
subroutine zero_densem(x,n,max_elem)
 ! Allocates and zero a dense matrix of n x n
 integer(i8) :: n
 integer(i8), optional :: max_elem !for compatibilty with zero sparse
 type (densem) :: x
!
 if (associated(x%x)) then
       if (n /= x%n) then
          deallocate(x%x)   ! different dimension, deallocate and allocate
          allocate(x%x(n,n))
       endif
    else
       allocate(x%x(n,n))   
 endif
 x%n=n
 x%x=0
end subroutine


subroutine zero_dense_symm(x,n,max_elem)
 ! Allocates and zero symmetric upper-stored dense matrix of n x n
 integer(i8) :: n
 integer(i8), optional :: max_elem !for compatibilty with zero sparse
 type (dense_symm) :: x
!
 if (associated(x%x)) then
       if (n /= x%n) then
          deallocate(x%x)   ! different dimension, deallocate and allocate
          allocate(x%x((n*(n+1))/2))
       endif
    else
       allocate(x%x((n*(n+1))/2))
 endif
 x%n=n
 x%x=0
end subroutine

subroutine zero_sparse_hashm(x,n,max_elem,storage_type)
 ! Allocates and zero sparse matrix in hash form
 integer(i8) :: n
 integer(i8), optional :: max_elem
 type (sparse_hashm) :: x
 character,optional::storage_type
!
!
 if (present(max_elem)) then
        x%nel=max_elem 
        x%nel=2** (1+ int(log(real(x%nel-1,r8))/log(2d0),i8)) ! AL that is, find the closest, higher, power of two
	                                                   ! 5 -> 8 ; 8 -> 8
							   ! http://stackoverflow.com/questions/4398711/round-to-the-nearest-power-of-two/4398995#4398995
	!print *,max_elem,x%nel,'--------'
    else
           ! for matrix used previously apply old dimensions
        if (.not. associated(x%xind)) x%nel=default_hash_size
 endif
!
 if (associated(x%xind)) then
       if (n /= x%n) then
          deallocate(x%xind)   ! different dimension, deallocate and allocate
          allocate(x%xind(2,x%nel))
          deallocate(x%xreal)   ! different dimension, deallocate and allocate
          allocate(x%xreal(1,x%nel))
       endif
    else
       allocate(x%xind(2,x%nel))
       allocate(x%xreal(1,x%nel))
 endif
 x%n=n; x%filled=0; x%status=1
 x%xreal=0
 x%xind=0
 x%storage_type='u'
 if (present(storage_type)) then
     if (storage_type == 'f') x%storage_type='f'
 endif
end subroutine

subroutine zero_sparse_ija(x,n)
 ! Allocates and zero sparse matrix in ija form
 integer(i8) :: n
 type (sparse_ija) :: x
!
! Actual allocation for the IJA format can be done through conversion
! from other programs, so here only deallocate if allocated

 if (associated(x%ia)) then
      call reset(x)
 endif
end subroutine


subroutine zero_ll(x,n)
 ! Allocates and zero a linked-list matrix of n x n
 integer(i8) :: n,i
 type (sparse_ll) :: x
 type(ll_element),pointer::current
!
 if (associated(x%ia)) then
       
       if (n /= x%n) then
            ! different dimension, deallocate and allocate
            ! call reset(x)
             print*,'reset in sparse_ll not yet implemented'
             stop
          else
           ! zero all elements but do not deallocate
             do i=1,n
                current=>x%ia(i)
                do while (associated(current))
                   current%val=0
                   current=>current%next
                enddo  
             enddo
       endif
    else
       ! new matrix
       allocate(x%ia(n))
       nullify(x%ia)
 endif
 x%n=n
end subroutine


!----------------------------
!  reset (deallocate) matrix
!----------------------------

subroutine deall_densem(x)
 ! Deallocates  dense matrix
 type (densem) :: x
!
 x%n=0
 if (associated(x%x)) deallocate(x%x)
end subroutine

subroutine deall_dense_symm(x)
 ! deallocates  symmetric upper-stored dense matrix 
 type (dense_symm) :: x
!
 x%n=0
 if (associated(x%x)) deallocate(x%x)
end subroutine

subroutine deall_sparse_hashm(x)
 ! deallocates sparse hash matrix
 type (sparse_hashm) :: x
!
 x%n=0; x%filled=0; x%status=0
 x%nel=0; x%storage_type='u'
 if (associated(x%xreal)) deallocate(x%xreal)
 if (associated(x%xind)) deallocate(x%xind)
end subroutine

subroutine deall_sparse_ija(x)
 ! deallocates sparse ija matrix
 type (sparse_ija) :: x
!
 x%n=0; x%nel=0
 if (associated(x%ia))deallocate(x%ia,x%ja,x%a)
end subroutine

!--------------
!  add to matrix
!--------------

subroutine add_densem(a,i,j,x)
! adds a(i,j) to dense matrix x
 real (rh) :: a
 type (densem):: x
 integer(i8) i,j

 x%x(i,j)=x%x(i,j)+a
end subroutine


subroutine add_dense_symm(a,i,j,x)
! adds a(i,j) to dense symmetric upper-stored matrix x
 real (rh) :: a
 type (dense_symm):: x
 integer(i8) i,j,k
 if (j >=i) then
    k=pos_symm_sp(i,j,x%n)
    x%x(k)=x%x(k)+a
 else
    stop 'Warning !!! this is an upper-triangular matrix'
 endif
end subroutine


 function pos_symm_sp(i,j,n) result (address)
 !finds position of (i,j) element of a n.n matrix upper-stored
 !renamed for compiler compatibility with an identical function in denseop module
 integer(i8) :: i,j,n,address
 
 if (j >= i) then
      address=(i-1)*n-(i*(i-3))/2+j-i
   else
      address=(j-1)*n-(j*(j-3))/2+i-j
 endif
 end function


recursive subroutine add_sparse_hashm(a,i,j,x)
! adds a(i,j) to sparse hash upper-stored matrix x
! storage_type is 'u' for upper-trangular (deafult) and 'f' for full
! get it from sparse_hashm structure 
 real (rh):: a
 type (sparse_hashm):: x,y
 integer(i8) i,j,k,full
 logical::upper_storage

 upper_storage=.true. 
 if (x%storage_type == 'f') upper_storage=.false.
     
 if (x%status /=1) then
    print*,'Structure sparse_hashm not ready; probably destroyed by solving'
    stop
 endif
 if (.not. upper_storage .or. j >=i) then
    !full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,1,x%filled)
    full=hashvr((/i,j/),2_i8,x%xind,x%nel,3_i8,1_i8,x%filled) !AL
    if(full < -1 ) print *,i,j,full
    if (full == -1 .or. real(x%filled,r8)/x%nel > maxfill ) then
!     Matrix full, copy to a matrix hash_incr times larger
       !print*,'x'
       !call printm(x,'internal')
       call init(y)
       !call zerom(y,x%n,int(hash_incr*x%nel,i8),x%storage_type)
       call zerom(y,x%n,2_i8*x%nel,x%storage_type)
       do k=1,x%nel
          if (x%xind(1,k) /= 0) then
              call addm(x%xreal(1,k),int(x%xind(1,k),i8),int(x%xind(2,k),i8),y)
          endif
       enddo
       !print*,'y'
       !call printm(y,'internal')
       print'(a,i0,a,i0,a,f10.4)','hash matrix increased from ',x%nel,' to ',2_i8*x%nel,&
                                  ' % filled: '!,real(x%filled,r8)/x%nel
       !call flush()                                  
       ! Move y to x
       call reset(x)
       x%n=y%n; x%nel=y%nel;x%filled=y%filled;x%status=y%status;x%xind=>y%xind; x%xreal=>y%xreal
       x%storage_type=y%storage_type
       nullify(y%xind)
       nullify(y%xreal)
       !print*,'new x'
       !call printm(x,'internal')
       !full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,1,x%filled)
       full=hashvr((/i,j/),2_i8,x%xind,x%nel,3_i8,1_i8,x%filled) !AL
    endif
    if (full >0) then
       x%xreal(1,full)=x%xreal(1,full)+a
       else
         call add_sparse_hashm(a,i,j,x)
    endif 
 endif
end subroutine

!--------------
!  set matrix element
!--------------

subroutine set_densem(a,i,j,x)
! x(i,j)=a for dense square matrix
 real (rh) :: a
 type (densem):: x
 integer(i8) i,j

 x%x(i,j)=a
end subroutine


subroutine set_dense_symm(a,i,j,x)
! x(i,j)=a for dense symmetric upper-stored matrix x
 real (rh) :: a
 type (dense_symm):: x
 integer(i8) i,j,k
 if (j >=i) then
    k=pos_symm_sp(i,j,x%n)
    x%x(k)=a
 else
    stop 'Warning !!! this is an upper-triangular matrix'
 endif
end subroutine


subroutine set_sparse_hashm(a,i,j,x)
! x(i,j)=a for sparse hash upper-stored or full-stored  matrix x

 real (rh) :: a
 type (sparse_hashm):: x
 integer(i8) i,j,full
 logical :: upper_storage
!
 upper_storage=.true.
 if (x%storage_type == 'f') upper_storage=.false.

 if (x%status /=1) then
    print*,'Structure sparse_hashm not ready; probably sorted'
    stop
 endif
 if (.not. upper_storage .or. j >=i) then
    !full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,1,x%filled)
    full=hashvr((/i,j/),2_i8,x%xind,x%nel,3_i8,1_i8,x%filled) !AL
    if (full == -1) then
       print*,'hash matrix too small'
       stop
    endif
    x%xreal(1,full)=a
 endif
end subroutine


!-------------------------
!  get scalar from matrix
!-------------------------

function get_densem(i,j,x) RESULT (a)
! a=x(i,j) for dense matrix x
 real (rh) :: a
 type (densem):: x
 integer(i8) i,j

 a=x%x(i,j)
end function


function get_dense_symm(i,j,x) RESULT(a)
! a=x(i,j) for dense symmetric upper-stored matrix x
 real (rh) :: a
 type (dense_symm):: x
 integer(i8) i,j,k
 if (j >=i) then
    k=pos_symm_sp(i,j,x%n)
    a=x%x(k)
 else
    stop 'Warning !!! this is an upper-triangular matrix'
 endif
end function


function get_sparse_hashm(i,j,x) RESULT(a)
! a=x(i,j) for sparse hash matrix

 real (rh) :: a
 type (sparse_hashm):: x
 integer(i8) i,j,full
 logical::upper_storage
!
 upper_storage=.true.
 if (x%storage_type == 'f') upper_storage=.false.

 if (x%status /=1) then
    print*,'Structure sparse_hashm not ready; probably sorted'
    stop
 endif
 if (.not. upper_storage .or. j >=i) then
    !full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,0,x%filled)
    full=hashvr((/i,j/),2_i8,x%xind,x%nel,3_i8,0_i8,x%filled) !AL
    if (full <1) then
       a=0
     else
       a=x%xreal(1,full)
    endif
  else
    a=0
 endif
end function

function get_sparse_ija(i,j,x) RESULT(a)
! a=x(i,j) for sparse IJA matrix
 real (r8) :: a
 type (sparse_ija) :: x
 INTEGER(i8),INTENT(IN) :: i,j
 integer(i8):: k

 do k=x%ia(i),x%ia(i+1)-1
   if (x%ja(k) ==j) then
       a=x%a(k)
       return
   end if
 enddo
 a=0    ! element not found
END function

end module
