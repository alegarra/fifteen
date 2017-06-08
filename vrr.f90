!VRR: Variance of Realized Relationships
! at one locus

! This programs computes Karigl's relationships and the covariances of realized
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
! modified to compute probabilities of being .not. IBD (i.e. heterogeneous) 6 jun 2017


program vrr
use kinds; use sparsem
implicit none
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !=real(8) con DVF
integer:: nanim,nsubset !=20 !22 !19 !14
integer:: nele
integer,allocatable:: ped(:,:),subset(:)
integer:: i,j,k,l,ii,jj,kk,t,io
integer::s1,s2,d1,d2
real(dp),allocatable:: F(:)
real(dp):: delta(15),val,bigF(4)
real(dp):: lhs(9,9),rhs(9),delta9(9),delta15(15),rhs15(15)
real(dp),allocatable:: vara(:,:),vara4d(:,:,:,:)
character(100):: pedfile,dum
logical:: use_hash(4)=.true.
integer(i16):: iab,iabc,iabpcd,iabcd
logical:: compute_var=.false.,compute_covar=.false.,compute_vardom=.false.,compute_dom=.false., &
                compute_iden=.false.,&
                compute_geninb=.false.,compute_inbcovar=.false.,colleau=.false.,MeuwLuo=.true.

type(sparse_hashm) :: phiab,phiabc,phiabcd,phiabpcd ! phi2, phi3, phi4, phi22


iab=0;iabc=0;iabpcd=0;iabcd=0
print *,'pedfile?'
read *,pedfile
print *,'use hash? T/F (4 answers needed for phi2 phi3 phi4 phi22)'
read *,use_hash
print *,'use Colleau? T/F '
read *,colleau
print *,'compute inbreeding and generalized inbreedings for some animals? T/F'
read(*,*)compute_geninb
print *,'compute regular inbreeding and their covariances for some animals?'
read(*,*)compute_inbcovar
print *,'compute coancestries and variances of coancestries ONLY of some animals? T/F'
read(*,*)compute_var
print *,'compute coancestries and variances AND covariances some animals? T/F'
read(*,*)compute_covar
!print *,'compute dominance coefficients of some animals? T/F'
print *,'compute dominance (d_R) coefficients of some animals? T/F'
read(*,*)compute_dom
print *,'compute variance of dominance coefficients of some animals? T/F'
read(*,*)compute_vardom
print *,'compute probabilities of identity coefficients? T/F'
read(*,*) compute_iden

!open(unit=1,file='lopetegui.txt',status='old')
!open(unit=1,file='mencha2.txt',status='old')
!open(unit=1,file='mencha.txt',status='old')
open(unit=1,file=pedfile,status='old')
!open(unit=4,file=adjustl(trim(pedfile))//'_results',status='replace')
nanim=0
do
        read(1,*,iostat=io)i
        if(io/=0) exit
        nanim=nanim+1
enddo        
rewind(1)
print *,'nanim= ',nanim; nele=nanim
allocate(F(0:nanim),ped(nanim,3))
call zerom(phiab,nanim,nele)
call zerom(phiabc,nanim,nele)
call zerom(phiabcd,nanim,nele)
call zerom(phiabpcd,nanim,nele)

do i=1,nanim
        read(1,*)ped(i,:)
enddo
close(1)

if(colleau) then
        F=0
        if(.not.MeuwLuo) then
                colleau=.false.
                do i=1,nanim
                        !if(mod(i,max((nanim/10),1))==0 ) print *,int(10*i/(nanim/10d0)),'%'
                        F(i)=2*phi2(i,i)-1
                enddo
                colleau=.true.
        else
                call meuw(ped(:,2),ped(:,3),f)
        endif
        print *,'inb done'
endif
if(compute_geninb) then
        call get_subset('generalized inbreeding')
        open(unit=4,file=adjustl(trim(dum))//'_geninb',status='replace')
        write(4,'(a)') 'i phi2 phi3 phi4 phi22'
        do i=1,nsubset
                ii=subset(i)
                if(mod(i,max((nsubset/10),1))==0 ) print *,int(10*i/(nsubset/10d0)),'%'
                bigF(1)=phi2(ii,ii)
                bigF(2)=phi3(ii,ii,ii)
                bigF(3)=phi4(ii,ii,ii,ii)
                bigF(4)=phi22(ii,ii,ii,ii)
                write(4,'(i10,1x,20f10.6)') ii,bigF
                
        enddo
        print *,'geninb done'
        close(4)
endif

if(compute_inbcovar) then
        call get_subset('regular inbreeding and covariances')
        open(unit=4,file=adjustl(trim(dum))//'_inbcovar',status='replace')
        write(4,'(a)') 'i j inb covar'
        do i=1,nsubset
                ii=subset(i)
                if(mod(i,max((nsubset/10),1))==0 ) print *,int(10*i/(nsubset/10d0)),'%'
                do j=i,nsubset
                        jj=subset(j)
                        if(ii==jj) then
                                val=2*phi2(ii,ii)-1
                        else
                                val=-99d0
                        endif
                        write(4,'(2i10,1x,20f10.6)') &
                          ii,jj,val,4*(phi22(ii,ii,jj,jj)-phi2(ii,ii)*phi2(jj,jj))
                enddo
                
        enddo
        print *,'inbcovar done'
        close(4)
endif
!this is the inverse of the matrix of coefficients in Karig (1981, eq.7)
lhs=reshape( (/ &
    0.,     0.,    0.,    1.,   - 1.,  - 1.,    1.,    0.,    0.,  &
    4.,   - 4.,  - 4.,  - 1.,     1.,    1.,  - 1.,    4.,    0.,  &
    0.,     0.,    0.,  - 4.,     4.,    2.,  - 2.,    0.,    0.,  & 
  - 8.,     8.,    4.,    4.,   - 4.,  - 2.,    2.,  - 4.,    0.,  &
    0.,     0.,    0.,  - 4.,     2.,    4.,  - 2.,    0.,    0.,  &
  - 8.,     4.,    8.,    4.,   - 2.,  - 4.,    2.,  - 4.,    0.,  &
    0.,     0.,    0.,    0.,     0.,    0.,  - 2.,    0.,    2.,  &
    0.,     0.,    0.,    16.,  - 8.,  - 8.,    8.,    0.,  - 4.,  &
    16.,  - 8.,  - 8.,  - 16.,    8.,    8.,  - 6.,    4.,    2. /),&
    (/9,9/) )
lhs=transpose(lhs)/4d0
! Miguel &A
! get identity coefficients 
if(compute_iden) then
        call get_subset('probabilities of identity states')
        open(unit=4,file=adjustl(trim(dum))//'_iden',status='replace')
        write(4,'(a)') 'i  j delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 delta9'
        do i=1,nsubset
                if(mod(i,max((nsubset/10),1))==0 ) print *,int(10*i/(nsubset/10d0)),'%'
                ii=subset(i)
                do j=i,nsubset
                        jj=subset(j)
                        rhs=(/ &
                                1d0,&
                                2*phi2(ii,ii),&
                                2*phi2(jj,jj),&
                                4*phi2(ii,jj),&
                                8*phi3(ii,ii,jj),&
                                8*phi3(ii,jj,jj),&
                                16*phi4(ii,ii,jj,jj),&
                                4*phi22(ii,ii,jj,jj),&
                                16*phi22(ii,jj,ii,jj) /)

                        Delta9=matmul(lhs,rhs)
                        write(4,'(2i10,1x,20f10.6)')ii,jj,delta9
                enddo
        enddo
        close(4)
        print *,'identities done'
endif

if(compute_dom)then
        ! ahora calculo de coeficientes de dominancia para le matriz
        ! de dominancia "no inbreeding" de DeBoer & Hoeschele 
        ! esta matriz tiene las probabilidades de estados de identidad delta9 y delta12, que sumqdos
        ! corresponden al estado condensado 7 tal como se calcula aqui
        ! es decir, los estados son 1212 y 1221
        call get_subset('dominance coefficients')
        open(unit=4,file=adjustl(trim(dum))//'_dom',status='replace')
        write(4,'(a)') 'i j dominance additive_relationship'
        do i=1,nsubset
                if(mod(i,max((nsubset/10),1))==0 ) print *,int(10*i/(nsubset/10d0)),'%'
                print *,i
                ii=subset(i)
                do j=1,nsubset
                    !if(mod(j,100)==0) print *,i,j
                    jj=subset(j)
                    val=dom(ii,jj)
                    !val=dom_R(ii,jj)
                    !if (val>0) then
                        ! only non-zero values are stored
                        ! individuals, dominance, additive relationship
                        if (i<=j) write(4,'(2i10,1x,3g20.12)')ii,jj,val,2d0*phi2(ii,jj)
                    !endif
                enddo
        enddo
        close(4)
        print *,'dominance finished'
endif

! vara4d
if(compute_covar)then
        call get_subset('covariance of coancestries')
        open(unit=4,file=adjustl(trim(dum))//'_covar',status='replace')
        write(4,'(a)') 'n nn i j k l cov(coan(i,j),coan(k,l)) coan(i,j) coan(k,l)'
        k=0
        do i=1,nsubset
                if(mod(i,max((nsubset/10),1))==0 ) print *,int(10*i/(nsubset/10d0)),'%'
                do ii=i,nsubset
                        k=k+1
                        kk=0
                        do j=1,nsubset
                                do jj=j,nsubset
                                        kk=kk+1
                                        !print *,subset(i),subset(ii),subset(j),subset(jj),iabpcd
                                        val=phi22(subset(i),subset(ii),subset(j),subset(jj))- &
                                                phi2(subset(i),subset(ii))*phi2(subset(j),subset(jj))
                                        write(4,'(6i10,1x,3g20.12)')k,kk,subset(i),subset(ii),subset(j),subset(jj),&
                                                        val,phi2(subset(i),subset(ii)),phi2(subset(j),subset(jj))
                                enddo
                        enddo
                enddo
        enddo
        print *,'covariances finished'
        close(4)
endif

if(compute_var)then
        call get_subset('variance of coancestries')
        open(unit=4,file=adjustl(trim(dum))//'_var',status='replace')
        write(4,'(a)') 'i j var(coan(i,j)) coan(i,j)'
        do i=1,nsubset
                if(mod(i,max((nsubset/10),1))==0 ) print *,int(10*i/(nsubset/10d0)),'%'
                do ii=i,nsubset
                                val=phi22(subset(i),subset(ii),subset(i),subset(ii))- &
                                        phi2(subset(i),subset(ii))*phi2(subset(i),subset(ii))
                                write(4,'(2i10,1x,3g20.12)')subset(i),subset(ii),&
                                                val,phi2(subset(i),subset(ii))
                enddo
        enddo
        print *,'variances finished'
        close(4)
endif

if(compute_vardom)then
        call get_subset('variance of dominances')
        open(unit=4,file=adjustl(trim(dum))//'_vardom',status='replace')
        write(4,'(a)') 'i j var_dominance dominance'
        do i=1,nsubset
                if(mod(i,max((nsubset/10),1))==0 ) print *,int(10*i/(nsubset/10d0)),'%'
                ii=subset(i)
                do j=i,nsubset
                        jj=subset(j)
                        rhs=(/ &
                                1d0,&
                                2*phi2(ii,ii),&
                                2*phi2(jj,jj),&
                                4*phi2(ii,jj),&
                                8*phi3(ii,ii,jj),&
                                8*phi3(ii,jj,jj),&
                                16*phi4(ii,ii,jj,jj),&
                                4*phi22(ii,ii,jj,jj),&
                                16*phi22(ii,jj,ii,jj) /)
                        Delta9=matmul(lhs,rhs)
                        val=Delta9(1)*(1-Delta9(1))+Delta9(7)*(1-Delta9(7))+2*Delta9(1)*Delta9(1)
                        write(4,'(2i10,1x,3g20.12)')ii,jj,&
                                                val,Delta9(1)+Delta9(7)
                enddo
        enddo
        print *,'variances of dominance finished'
        close(4)
endif

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

function phi2_colleau(k,l) result(s)
    implicit none
    integer(i8):: i,k,l
    real(r8):: s,w(nanim),v(nanim)
    v=0d0; v(k)=1d0
    call A_times_v(w,v)
    s=.5*w(l) !the .( is to pass from additive relationships to coancestries
end function    
    
subroutine A_times_v(w,v)
! computes w=A*v using  TDT'v
integer :: i
real(r8) :: w(:),v(:),tmp,di
real(r8) ::q(size(w))
q=0

do i=nanim,1,-1
   q(i)=q(i)+v(i)
   if (ped(i,2)>0) q(ped(i,2)) = q(ped(i,2))+q(i)*0.5
   if (ped(i,3)>0) q(ped(i,3)) = q(ped(i,3))+q(i)*0.5
   !q(ped(i,:)) = q(ped(i,:))+q(i)*0.5
enddo
do i=1,nanim
   di=(count(ped(i,2:3)==0)+2)/4d0 - .25*(f(ped(i,2))+f(ped(i,3)))
   tmp = 0
   if (ped(i,2)>0) tmp = tmp +w(ped(i,2))
   if (ped(i,3)>0) tmp = tmp +w(ped(i,3))
   w(i)=.5*tmp
   w(i)= w(i) + di*q(i)
enddo
end subroutine


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
                        s=0.2d0 ! this is actually gamma/2 !!!
                else
                    if(colleau) then
                        s=phi2_colleau(a,b)
                    else
                        if(a==b) then
                                s=.5*(1+phi2(ped(a,2),ped(a,3)))
                        else
                                s=.5*(phi2(ped(a,2),b)+phi2(ped(a,3),b))
                        endif
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
        pos1=(a+1)*nanim+b+1
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
                        s=0d0
                ! all equal aaa
                else if(all(a==(/b,c/))) then
                        s=.25*(1+3*phi2(ped(a,2),ped(a,3)))
                !aac
                else if(a==b) then
                        s=.5*(phi2(a,c)+phi3(ped(a,2),ped(a,3),c))
                else
                !abc
                        s=.5*(phi3(ped(a,2),b,c)+phi3(ped(a,3),b,c))
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
        pos1=(a+1)*nanim+b+1
        pos2=(c+1)*nanim+d+1
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
                        s=0d0
                ! all equal aaaa
                else if(all(a==(/b,c,d/))) then
                        s=.125*(1+7*phi2(ped(a,2),ped(a,3)))
                !aaad
                else if(all(a==(/b,c/)  )) then
                        s=.25*(phi2(a,d)+3*phi3(ped(a,2),ped(a,3),d))
                !aacd
                else if(all(a==(/b/)    )) then
                        s=.5*(phi3(a,c,d)+phi4(ped(a,2),ped(a,3),c,d))
                else
                !abcd
                        s=.5*(phi4(ped(a,2),b,c,d)+phi4(ped(a,3),b,c,d))
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
        pos1=(a+1)*nanim+b+1
        pos2=(c+1)*nanim+d+1
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
                s=0d0
        else 
                ! all equal aa,aa
                if(all(a==(/b,c,d/))) then
                        s=.25*(1+3*phi2(ped(a,2),ped(a,3)))
                else
                        if(a<b) call swap(a,b)
                        if(c<d) call swap(c,d)
                        if((a==c).and.(b<d)) call swap(b,d)
                        if(all(a==(/b,c/)  )) then
                                s=.5*(phi2(a,d)+phi3(ped(a,2),ped(a,3),d))
                        else
                                if(a<c) then
                                        call swap(a,c); call swap(b,d)
                                endif
                                if(a==b) then
                                        s=.5*(phi2(c,d)+phi22(ped(a,2),ped(a,3),c,d))
                                else
                                        if(a==c) then
                                                s=.25*(2*phi3(a,b,d)+phi22(ped(a,2),b,ped(a,3),d) &
                                                        +phi22(ped(a,3),b,ped(a,2),d))
                                        else
                                                s=.5*(phi22(ped(a,2),b,c,d)+phi22(ped(a,3),b,c,d))
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
        pos1=(a+1)*nanim+b+1
        pos2=(c+1)*nanim+d+1
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


function dom(a,b) result(s)
        implicit none
        integer(i8)::a,b
        real(dp):: s,rhs(9),Delta9(9)
        !lhs is in the main program
        rhs=(/ &
                1d0,&
                2*phi2(a,a),&
                2*phi2(b,b),&
                4*phi2(a,b),&
                8*phi3(a,a,b),&
                8*phi3(a,b,b),&
                16*phi4(a,a,b,b),&
                4*phi22(a,a,b,b),&
                16*phi22(a,b,a,b) /)
        Delta9=matmul(lhs,rhs)
        s=Delta9(1)+Delta9(7)
end function

function dom_R(a,b) result(s)
        implicit none
        integer(i8)::a,b
        real(dp):: s,rhs(9),Delta9(9)
        !lhs is in the main program
        rhs=(/ &
                1d0,&
                2*phi2(a,a),&
                2*phi2(b,b),&
                4*phi2(a,b),&
                8*phi3(a,a,b),&
                8*phi3(a,b,b),&
                16*phi4(a,a,b,b),&
                4*phi22(a,a,b,b),&
                16*phi22(a,b,a,b) /)
        Delta9=matmul(lhs,rhs)
        s=Delta9(7)
end function


subroutine meuw(s1,d1,f) 
   !  Meuwissen & Luo
   implicit none
   integer(i8),intent(in):: s1(:),d1(:)
   integer(i8) :: ss(size(s1)),dd(size(s1))
   real(dp):: f(0:)
   integer(i8), allocatable:: point(:)
   real(r8), allocatable:: l(:),d(:)
   integer(i8) :: is,id,i,j,k,n,ks,kd
   real(r8) :: r,fi
   real :: t1

   print*,'Calculating Inbreeding by M&L function'
   t1=seconds()
   ss=s1; dd=d1
   n=size(ss)
   allocate(point(n),l(n),d(n))
   point=0;l=0;d=0
   f(0)=-1
   do i=1,n
      if(n>10 .and. mod(i,int(n/10))==0) &
        write(*,*) ' at',i,' animals'
      is=ss(i); id=dd(i)
      ss(i)=max(is,id)
      dd(i)=min(is,id)
      d(i)=.5-.25*(f(is)+f(id))
      if (is==0 .or. id==0) then
         f(i)=0
         cycle
      endif
      fi=-1
      l(i)=1.0
      j=i
      do while (j/=0)
         k=j
         r=0.5*l(k)
         ks=ss(k)
         kd=dd(k)
         if (ks>0) then
            do while(point(k)>ks)
               k=point(k)
            enddo
            l(ks)=l(ks)+r
            if (ks/=point(k))then
               point(ks)=point(k)
               point(k)=ks
            endif
            if (kd>0)then
               do while (point(k)>kd)
                  k=point(k)
               enddo
               l(kd)=l(kd)+r
               if (kd/=point(k)) then
                  point(kd)=point(k)
                  point(k)=kd
               endif
            endif
         endif

         fi=fi+l(j)*l(j)*d(j)
         l(j)=0
         k=j
         j=point(j)
         point(k)=0
      enddo
      f(i)=fi
   enddo
   f(0)=0
   print*,'   Calculating Inbreeding by M&L, elapsed time',seconds()-t1
end subroutine



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

end


