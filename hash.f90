module hash
use kinds
public::mix,final,rot,hashvr
contains
!Fortran implementation of lookup3.c from Bob Jenkins
!http://burtleburtle.net/bob/hash/index.html#lookup
!Converted in fortran by Francois Guillaume- 2011
  
integer(kind=i8) function rot(i,j)
  !======================================================================
  !Rotation function
  !======================================================================
  use kinds
  implicit none
  integer(kind=i8) :: i,j
  rot=IOR(ISHFT(i,j),ISHFT(i,-(32_i8-j)))
  return
end function rot
    
subroutine mix(a,b,c)
  !======================================================================
  !Mixing subroutine
  !======================================================================
  use kinds
  implicit none
  integer(kind=i8),intent(inout)::a,b,c
  a=a-c  ; a=IEOR(a,rot(c, 4_i8)) ;c=c+b 
  b=b-a  ; b=IEOR(b,rot(a, 6_i8)) ;a=a+c 
  c=c-b  ; c=IEOR(c,rot(b, 8_i8)) ;b=b+a 
  a=a-c  ; a=IEOR(a,rot(c,16_i8)) ;c=c+b 
  b=b-a  ; b=IEOR(b,rot(a,19_i8)) ;a=a+c 
  c=c-b  ; c=IEOR(c,rot(b, 4_i8)) ;b=b+a 
end subroutine mix

subroutine final(a,b,c)
  !======================================================================
  !Final Mixing subroutine (for 4+ coordinates)
  !======================================================================  
  use kinds
  implicit none      
  integer(kind=i8),intent(inout)::a,b,c
  c=IEOR(c,b); c=c-rot(b,14_i8)
  a=IEOR(a,c); a=a-rot(c,11_i8)
  b=IEOR(b,a); b=b-rot(a,25_i8)
  c=IEOR(c,b); c=c-rot(b,16_i8)
  a=IEOR(a,c); a=a-rot(c,4_i8) 
  b=IEOR(b,a); b=b-rot(a,14_i8)
  c=IEOR(c,b); c=c-rot(b,24_i8)
end subroutine final


integer(i8) function hashvr(dat,r,ind,m,n,mode,nr)
  !======================================================================
  ! Hash function
  !======================================================================  
  ! returns address of  element of ind()  containing integer dat.
  !c Searching goes with a HASH technique. When mode=1 and dat was not
  !c there before, dat is written  into ind("hash",1..r). When mod=0
  !c and dat is not found, this function returns 0.
  !c rank of ind is m x n,
  !c dat contains r numbers. 
  !======================================================================  
  use kinds
  implicit none
  integer(i8) r,mode,nr,m,n
  integer(i8)  :: ind(:,:)
  !real(rh)  :: ind(:,:) !AL
  integer(i8) :: iaddr,k,izer,ieq,i1,i,dat(r)
  integer(kind=i8) ::a,b,c,iaddress,plage
  integer(i8) ::Mult,Reste,Passe
  
  !How many time should we mix coordinates
  Mult=int(r/3_i8);Reste=mod(r,3_i8)
  !Init
  Passe=0_i8;                     !How 
  plage=int(m,kind(a))         !size of the array
  a=int(dat(1),kind(a))        !Conversion of first coordinate
  b=-1640531527_i8             !Default value for 2nd coordinate 
  if(r>1)b=int(dat(2),kind(a)) !
  c=305419896_i8               !Default value for 3rd coordinate  
  if(r>2)c=int(dat(3),kind(a))
  !Cycling in order to hash all r values
  HashInit:do 
     call mix(a,b,c)
     if(Passe<Mult)then 
        Passe=Passe+1_i8   
        if(Passe<Mult)then 
           a=int(dat(1_i8+(3*Passe)),kind(a))  ! Un batch de 3
           b=int(dat(2_i8+(3*Passe)),kind(a)) ! entiers a faire 
           c=int(dat(3_i8+(3*Passe)),kind(a)) ! passer
           cycle HashInit
        else !on a un modulo non nul et derniere passe
           if(Reste/=0)then
              a=int(dat(1_i8+(3*Passe)),kind(a))
              if(Reste==2_i8)b=int(dat(2_i8+(3*Passe)),kind(a))
              call final(a,b,c)
              exit HashInit
           endif
        endif
     endif
     exit HashInit
  end do HashInit
  
  !Computation of the address
  iaddress=IAND(c, plage - 1_i8) + 1_i8
  if(iaddress<0_i8)then
        print *,dat,iaddress
        stop
  endif
  
  !Cycling until a free cell is found
  Hash:do  k=1,5000
     izer=0;ieq=0
     do i=1,r
        i1=ind(i,iaddress)
        if (i1.ne.dat(i)) ieq=1
        if (ind(i,iaddress).ne.0) izer=1
     enddo     
     if (izer.eq.0 .or. ieq.eq.0) then
        if (izer.eq.0 .and. mode.eq.1) then
           do  i=1,r
              ind(i,iaddress)=dat(i)
           enddo
           nr=nr+1
        endif
        if (mode.eq.0 .and.izer.eq.0) then
           hashvr=0
        else
           hashvr=iaddress
        endif
        if(hashvr<-1) print *,'hashvr',iaddress,'dat',dat(1:2)
        return
     endif
     !Hashing again
     call mix(a,b,c)
     iaddress=IAND(c, plage - 1_i8) + 1_i8
  enddo Hash
  hashvr=-1
  write(*,*)"Pb 5000 collisions"
  return 
end function hashvr
end module hash

