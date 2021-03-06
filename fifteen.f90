
! This programs computes LAGC gametic relationships
! Karigl's relationships and the covariances of realized
! relationships at one locus following our derivations (LAGC, MATI, CCh & AL)
! a hash storage allows not recomputing coefficients again and again
! yet it is memory consuming.
! Also can compute dominance coefficients using dom() 
! covariance of dominance coefficients not implemented

! extremely slow for Lacaune ped. A way of speeding it is to consider that 
! Finally, PHIabc = 0 (or PHI( a b c d) = 0) if there is no common ancestor of a , b ,
! c (and d ) and
! @ a b , c d = 0 in the absence of two common ancestors, one for a and b and one
! for c and d . 
! so, a way to compute phiabc, phi ab,cd would be to check if those ancestors
! exist or not. Not done

! the pedigree has to be renumbered. Parents need to be
! renumbered before offspring, and it has to be ordered 1...n with no holes
! missing parent is indicated as 0

! does a few tasks:
!       o compute variances (_var) and covariances of coancestries (file  _covar)
!       o compute variances of dominance coefficients ( _vardom)
!       o compute dominance coefficients                ( _dom)
!       o compute probabilities of identity states      ( _iden)
!       o compute generalized inbreeding coefficients  ( _geninb)

! all for a subset of animals a,d their combinations defined in another file
! A Legarra, 7/8/2012, get/addmAsp functions from I Aguilar, hash functions from blupf90 & F Guillaume, 
! phi functions from Karigl 1981
! psi functions from Garcia-Cortes 2015
! modified to compute probabilities of being .not. IBD (i.e. heterogeneous) 6 jun 2017

! TO DO: think how to include metafounders so far this is done including implictely 0.1 !!!


program vrr
use kinds; use sparsem
implicit none
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !=real(8) con DVF
integer:: nanim,nsubset !=20 !22 !19 !14
integer:: nele,nmf,ngam
integer,allocatable:: ped(:,:),subset(:),pedgam(:,:)
integer:: i,j,k,l,ii,jj,kk,t,io,sire,dam,i1,i2,j1,j2
integer::s1,s2,d1,d2
real(dp),allocatable:: F(:),gamma(:,:)
real(dp):: delta(15),val,bigF(4)
real(dp):: lhs(15,15),rhs(15),delta9(15),delta15(15),rhs15(15)
real(dp),allocatable:: vara(:,:),vara4d(:,:,:,:)
character(100):: pedfile,dum
logical:: use_hash(4)=.true.
integer(i16):: iab,iabc,iabpcd,iabcd

!type(sparse_hashm) :: phiab,phiabc,phiabcd,phiabpcd ! phi2, phi3, phi4, phi22
type(sparse_hashm) :: phiab,phiabc,phiabcd,phiabpcd ! phi2, phi3, phi4, phi22
! REPLACE PHI BY PSI
! rule: psi(a,a)=1; psi(0,0)=gamma

iab=0;iabc=0;iabpcd=0;iabcd=0
print *,'pedfile?'
read *,pedfile

! use hash, output 15 states only


lhs=reshape( (/ &
  1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
 -1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
 -1.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
 -1.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
 -1.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
 -1.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
  2.0,  -1.0,  -1.0,   0.0,   0.0,  -1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
  2.0,   0.0,   0.0,  -1.0,  -1.0,  -1.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
 -1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
  2.0,  -1.0,   0.0,  -1.0,   0.0,   0.0,   0.0,   0.0,  -1.0,   1.0,   0.0,   0.0,   0.0,   0.0,  0.0, &
  2.0,   0.0,  -1.0,   0.0,  -1.0,   0.0,   0.0,   0.0,  -1.0,   0.0,   1.0,   0.0,   0.0,   0.0,  0.0, &
 -1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,  0.0, &
  2.0,   0.0,  -1.0,  -1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  -1.0,   1.0,   0.0,  0.0, &
  2.0,  -1.0,   0.0,   0.0,  -1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  -1.0,   0.0,   1.0,  0.0, &
 -6.0,   2.0,   2.0,   2.0,   2.0,   1.0,  -1.0,  -1.0,   1.0,  -1.0,  -1.0,   1.0,  -1.0,  -1.0,  1.0 /),&
(/15,15/) )
lhs=transpose(lhs)

do i=1,15
    write(*,'(15f6.3)')lhs(i,:)
enddo

!open(unit=1,file='lopetegui.txt',status='old')
!open(unit=1,file='mencha2.txt',status='old')
!open(unit=1,file='mencha.txt',status='old')
open(unit=1,file=pedfile,status='old')
!open(unit=4,file=adjustl(trim(pedfile))//'_results',status='replace')
nanim=0
ngam=nanim*2
nmf=0
do
        read(1,*,iostat=io)i,j,k
        if(io/=0) exit
        nanim=nanim+1
        if(j==0) nmf=nmf+1
enddo        
ngam=nanim*2
rewind(1)
print *,'nanim= ',nanim
print *,'number of gametes= ',ngam
print *,'nmf= ',nmf
nele=ngam
allocate(F(nanim),ped(nanim,3),pedgam(ngam,3))
call zerom(phiab,ngam,nele)
call zerom(phiabc,ngam,nele)
call zerom(phiabcd,ngam,nele)
call zerom(phiabpcd,ngam,nele)

do i=1,nanim
        read(1,*)ped(i,:)
enddo
rewind(1)
!do i=1,nmf
!    read(1,*)ped(i,:),gamma(i,:)
!enddo
close(1)

! create gametic pedigree from animal pedigree
pedgam=0
k=0
do i=1,nanim
        ! paternal gamete
        k=k+1
        pedgam(k,1)=k
        ! flag sire
        sire=ped(i,2)
        ! expand sire to two gametes corresponding to the sire
        pedgam(k,2)=anim2gam(sire,1)
        pedgam(k,3)=anim2gam(sire,2)
        write (*,'(3i6,a,2i6,a,2i6,a,2i6,a,3i6)') pedgam(k,:),'|',gam2anim(pedgam(k,1)),';',gam2anim(pedgam(k,2)),';',gam2anim(pedgam(k,3)),'|',ped(i,:)
        ! maternal gamete
        k=k+1
        pedgam(k,1)=k
        ! flag dam
        dam=ped(i,3)
        ! expand sire to two gametes corresponding to the sire
        pedgam(k,2)=anim2gam(dam,1)
        pedgam(k,3)=anim2gam(dam,2)
        write (*,'(3i6,a,2i6,a,2i6,a,2i6,a,3i6)') pedgam(k,:),'|',gam2anim(pedgam(k,1)),';',gam2anim(pedgam(k,2)),';',gam2anim(pedgam(k,3)),'|',ped(i,:)
enddo

open(2,file='Dh.out',status='replace')
write(2,*) 'i j D'
do i=1,nanim
    do j=i,nanim
        !i1,i2,j1,j2
        i1=anim2gam(i,1)
        i2=anim2gam(i,2)
        j1=anim2gam(j,1)
        j2=anim2gam(j,2)
        rhs=0
        rhs(1)= phi4(i1,i2,j1,j2)
        rhs(2)= phi3(i1,i2,j1)
        rhs(3)= phi3(i1,i2,   j2)
        rhs(4)= phi3(i1,   j1,j2)
        rhs(5)= phi3(   i2,j1,j2)
        rhs(6)=phi22(i1,i2,j1,j2)
        rhs(7)= phi2(i1,i2)
        rhs(8)= phi2(j1,j2)
        rhs(9)=phi22(i1,j1,i2,j2)
        rhs(10)=phi2(i1,j1)
        rhs(11)=phi2(i2,j2)
        rhs(12)=phi22(i1,j2,i2,j1)
        rhs(13)=phi2(i1,j2)
        rhs(14)=phi2(i2,j1)
        rhs(15)=1
        Delta15=matmul(lhs,rhs)
        !write(*,'(2i4,20f9.5)') i,j,rhs
        write(*,'(2i4,20f9.5)') i,j,Delta15,sum(Delta15)
        write(2,*)i,j,1+phi22(i1,i2,j1,j2)-phi2(i1,i2)-phi2(j1,j2)!(the same),sum(Delta15(9:15))
    enddo
enddo



! information
print *,'total number of non-0 elements in phiab',phiab%filled
print *,'total number of non-0 elements in phiabc',phiabc%filled
print *,'total number of non-0 elements in phiabcd',phiabcd%filled
print *,'total number of non-0 elements in phiabpcd',phiabpcd%filled
print *,'total number of calls to phi2',iab
print *,'total number of calls to phi3',iabc
print *,'total number of calls to phi4',iabcd
print *,'total number of calls to phi22',iabpcd


contains


recursive double precision function phi2(k,l) result(s)
        implicit none
        integer(i8):: k,l,temp(2)
        integer(i8):: a,b
!        double precision:: s
        real(rh):: val
        integer(i8):: pos1,pos2
        iab=iab+1
        ! a must be younger (higher number) than b
        temp=bsort((/k,l/))
        a=temp(1); b=temp(2)
        ! already stored?
        pos1=a+1; pos2=b+1 ! the one is needed to make room for 0
        call getmAsp(val,pos1,pos2,phiab,use_hash(1))
        s=val
        ! how to differentiate an empty location from true:
                !empty = 0
                !0 is stored as -10
        if(s /= 0d0) then
                if(s == -10d0) s=0d0
                !print *,'ab'
                !return ! leave the subroutine
        else
                ! if it has not, compute s
                !if(a*b==0) then
                !        s=0d0
                if(a==0 .and. b==0) then
                        s=0.1d0 ! this is actually gamma/2 !!!
                else
                        if(a==b) then
                                s=1.d0
                                !s=.5*(1+phi2(pedgam(a,2),pedgam(a,3)))
                        else
                                s=.5*(phi2(pedgam(a,2),b)+phi2(pedgam(a,3),b))
                        endif
                endif
                val=s
                if(s==0d0) val=-10d0
                call addmAsp(real(val,rh),pos1,pos2,phiab,use_hash(1))
        endif
end function



recursive double precision function phi3(j,k,l) result(s)
        implicit none
        integer(i8),intent(in):: j,k,l
        integer(i8):: temp(3),a,b,c
        integer(i8):: pos1,pos2
        real(rh)::val
        iabc=iabc+1
        !if(mod(iabc,100000000_i8)==0) print *,'phi3 iabc',iabc
        !sort so j>k>l (a is older than the other two)
        temp=bsort( (/j,k,l/) )
        !print *,'temp',temp
        a=temp(1); b=temp(2); c=temp(3)
        ! already stored?
        pos1=(a+1)*ngam+b+1
        pos2=c+1
        call getmAsp(val,pos1,pos2,phiabc,use_hash(2))
        s=val
        ! how to differentiate an empty location from true:
                !empty = 0
                !0 is stored as -10
        if(s /= 0d0) then
                !print *,'abc'
                if(s == -10d0) s=0d0
                !return ! leave the subroutine
        else
                if(a*b*c==0) then
                        s=0.1d0
                ! all equal aaa
                else if(all(a==(/b,c/))) then
                        s=1d0
                        !s=.25*(1+3*phi2(pedgam(a,2),pedgam(a,3)))
                !aac
                else if(a==b) then
                        s=phi2(a,c)
                        !s=.5*(phi2(a,c)+phi3(pedgam(a,2),pedgam(a,3),c))
                else
                !abc
                        s=.5*(phi3(pedgam(a,2),b,c)+phi3(pedgam(a,3),b,c))
                endif
                val=s
                if(s==0d0) val=-10d0
                call addmAsp(real(val,rh),pos1,pos2,phiabc,use_hash(2))
        endif
end function

recursive double precision function phi4(j,k,l,m) result(s)
        implicit none
        integer(i8),intent(in):: j,k,l,m
        integer(i8):: temp(4),a,b,c,d
        integer(i8):: pos1,pos2
        real(rh)::val
        iabcd=iabcd+1
        !sort so j>k>l>m
        temp=bsort( (/j,k,l,m/) )
        !print *,'temp',temp
        a=temp(1); b=temp(2); c=temp(3); d=temp(4)
        ! already stored?
        pos1=(a+1)*ngam+b+1
        pos2=(c+1)*ngam+d+1
        call getmAsp(val,pos1,pos2,phiabcd,use_hash(3))
        s=val
        ! how to differentiate an empty location from true:
                !empty = 0
                !0 is stored as -10
        if(s /= 0d0) then
                !print *,'abcd'
                if(s == -10d0) s=0d0
                !return ! leave the subroutine
        else
                if(a*b*c*d==0) then
                        s=0.1d0
                ! all equal aaaa
                else if(all(a==(/b,c,d/))) then
                        s=1d0
                        !s=.125*(1+7*phi2(pedgam(a,2),pedgam(a,3)))
                !aaad
                else if(all(a==(/b,c/)  )) then
                        s=phi2(a,d)
                        !s=.25*(phi2(a,d)+3*phi3(pedgam(a,2),pedgam(a,3),d))
                !aacd
                else if(all(a==(/b/)    )) then
                        s=phi3(a,c,d)
                        !s=.5*(phi3(a,c,d)+phi4(pedgam(a,2),pedgam(a,3),c,d))
                else
                !abcd
                        s=.5*(phi4(pedgam(a,2),b,c,d)+phi4(pedgam(a,3),b,c,d))
                endif
                val=s
                if(s==0d0) val=-10d0
                call addmAsp(real(val,rh),pos1,pos2,phiabcd,use_hash(3))
        endif
end function

recursive double precision function phi22(j,k,l,m) result(s)
        implicit none
        integer(i8),intent(in):: j,k,l,m
        integer(i8):: temp(4),a,b,c,d
        integer(i8):: pos1,pos2
        !double precision:: s
        real(rh)::val
        ! I am following Karigl here because a,b,c,d is very MESSI ;-)
        iabpcd=iabpcd+1
        if(mod(iabpcd,100000000_i8)==0) print *,'phi22 iabpcd',iabpcd
        ! do a preliminary sort to check if the value is stored
        a=j; b=k; c=l; d=m
        if(a<b) call swap(a,b)
        if(c<d) call swap(c,d)
        pos1=(a+1)*ngam+b+1
        pos2=(c+1)*ngam+d+1
        call getmAsp(val,pos1,pos2,phiabpcd,use_hash(4))
        s=val
        ! how to differentiate an empty location from true:
                !empty = 0
                !0 is stored as -10
        if(s /= 0d0) then
                !print *,'abpcd'
                if(s == -10d0) s=0d0
                return ! leave the subroutine
        endif

        ! if it is absent
        a=j; b=k; c=l; d=m
        if((a*b==0).and.(c*d==0)) then
                s=0.1d0
        else 
                ! all equal aa,aa
                if(all(a==(/b,c,d/))) then
                        s=1d0
                        !s=.25*(1+3*phi2(pedgam(a,2),pedgam(a,3)))
                else
                        if(a<b) call swap(a,b)
                        if(c<d) call swap(c,d)
                        if((a==c).and.(b<d)) call swap(b,d)
                        if(all(a==(/b,c/)  )) then
                                s=phi2(a,d)
                                !s=.5*(phi2(a,d)+phi3(pedgam(a,2),pedgam(a,3),d))
                        else
                                if(a<c) then
                                        call swap(a,c); call swap(b,d)
                                endif
                                if(a==b) then
                                        s=phi2(c,d)
                                        !s=.5*(phi2(c,d)+phi22(pedgam(a,2),pedgam(a,3),c,d))
                                else
                                        if(a==c) then
                                                s=phi3(a,b,d)
                                                !s=.25*(2*phi3(a,b,d)+phi22(pedgam(a,2),b,pedgam(a,3),d) &
                                                !        +phi22(pedgam(a,3),b,pedgam(a,2),d))
                                        else
                                                s=.5*(phi22(pedgam(a,2),b,c,d)+phi22(pedgam(a,3),b,c,d))
                                        endif
                                endif
                        endif
                endif
        endif
        val=s
        if(s==0d0) val=-10d0
        ! get adresses back after so many switches
        a=j; b=k; c=l; d=m
        if(a<b) call swap(a,b)
        if(c<d) call swap(c,d)
        pos1=(a+1)*ngam+b+1
        pos2=(c+1)*ngam+d+1
        call addmAsp(val,pos1,pos2,phiabpcd,use_hash(4))
end function

subroutine swap(a,b)
implicit none
integer(i8)::a,b,c
c=a
a=b
b=c
end subroutine

subroutine addmAsp(ai,ii,jj,x,use_hash)
! does not work with 0 or negative coordinates !!
! a=x(i,j) for sparse hash matrix
! with internal swap of variables
   implicit none
   real(rh) ,intent(in) :: ai
   real(rh) :: a
   type (sparse_hashm):: x
   integer(i8) , intent(in) :: ii,jj
   integer(i8) ::swap,i,j
   logical:: use_hash
   !
   if(.not.use_hash) return
   a=ai
   i=ii
   j=jj
   if (ai==0.e0)  return
   if (i>j) then
      swap=i
      i=j
      j=swap
   endif
   call addm(a,i,j,x)
end subroutine

subroutine getmAsp(ai,ii,jj,x,use_hash)
! does not work with 0 or negative coordinates !!
! a=x(i,j) for sparse hash matrix
   implicit none
   real(rh) ,intent(out) :: ai
   real (rh) :: a
   integer(i8) , intent(in) :: ii,jj
   type (sparse_hashm):: x
   integer(i8) i,j,swap
   logical:: use_hash
   !
   if(.not.use_hash) then
        ai=0
        return
   endif
   i=ii
   j=jj
   if (i>j) then
      swap=i
      i=j
      j=swap
   endif
   a=getm(i,j,x)
   ai=a
end subroutine


function bsort(a) result(b)
! _____________________
! sort in reverse order !!!
! ---------------------
implicit none
integer(i8):: a(:),b(size(a)),x,i,n
logical:: changed

n=size(a)
b=a
do 
        changed=.false.
        do i=1,n-1
                if(b(i)<b(i+1)) then
                        call swap(b(i),b(i+1))
                        changed=.true.
                endif
        enddo
        if(.not.changed) exit
enddo
end function






subroutine get_subset(message)
implicit none
character(len=*):: message
        print *, 'file with subset of animals to compute ',trim(adjustl(message)), ' (''all'' means all)'
        read *,dum
        if (allocated(subset)) deallocate (subset)
        if (dum=='all') then
                dum=pedfile
                nsubset=nanim
                allocate(subset(nsubset))
                subset=ped(:,1)
        else
                open(1,file=dum,status='old')
                nsubset=0
                do
                        read(1,*,iostat=io) i
                        if(io/=0) exit
                        nsubset=nsubset+1
                enddo
                print *,'size of the subset: ',nsubset
                allocate(subset(nsubset))
                rewind(1)
                do i=1,nsubset
                        read(1,*,iostat=io) subset(i)
                enddo
                close(1)
        endif
                
end subroutine

 real function seconds()
     call cpu_time(seconds)
      end function

integer function anim2gam(animal,origin) !origin is 1,2 (paternal, maternal)
    implicit none
    integer::animal,origin
    anim2gam=(animal-1)*2+origin
    if (anim2gam<0) anim2gam=0
    end function

 function gam2anim(gamete) result(pair) ! on output animal, pair
 implicit none
 integer::gamete,pair(2)
 pair(1)=gamete/2+mod(gamete,2) 
 pair(2)=gamete-2*(pair(1)-1)
 end function

end




