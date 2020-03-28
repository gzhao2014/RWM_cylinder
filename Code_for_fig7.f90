MODULE DATA_MODULE ! 参数模块 For kai eigenvalue
    INTEGER(4),PARAMETER:: N=2,miu=1,nn=500 
                                                               !number of equations
    COMPLEX(8),PARAMETER:: j=(0.,1.)      !pressure in injoined using integration method                                                                ! REAL(8),PARAMETER:: v0,b,k,m,d,deta, 
    REAL(8),PARAMETER::PI=3.1415926_8,r0=10,bb0=1.8,p0=4.0e5,miu0=12.56e-7,er=1.0e-20,tol=0.1,hmax=1.0e-3,hmin=1.0e-5 !B  hmin for the derivation of coefficients 
    END MODULE DATA_MODULE !p0=4.8 4.0e4aa=0.66*2/10

   PROGRAM gen      !PRINT THE value directly and scan the qa value for different parameters 
    
    USE DATA_MODULE    
	IMPLICIT NONE                                                                                 ! USE BESSI1，BESSK1，BESSI，BESSK !solve the bessel functiong 1 order and higher order
   real(8):: check,wou,wii,wii0,a,w2,w1,vc,aa,da,v0,v1,k,k11,d,bsita,bt,h,b,yita,yita1,dakai,dakai1,a1,rr,deta,bessi,eps,dw,apha,bb3,qa,r,q0,qq,wiu,wru,wr,bb3i,dwr,dv                        !dd=denominator of c1... 
                                                                                                       !h--step length, deta --> conductive rate,d -->shell thickness den--density bt toroidal field
                                                                                                       !gder=derivative of g,
   COMPLEX(8),DIMENSION(N):: Z,K1,K2,K3,K4,k5,k6,zout
   INTEGER(4):: i,m,jj,kj,jk,jkn,bn      
   COMPLEX(8),DIMENSION(nn):: wv
   COMPLEX(8):: w,mm,alph,ws,beta
   COMPLEX(8):: bb1,bb2,bb2u,z2h,z22h,wei,weider,z1h,wei0,wei00
   external  bessi
      !OPEN(UNIT=10,file="moide5.dat")
      OPEN(UNIT=1,file="eigenvalue.dat")
      !OPEN(UNIT=2,file="wr.dat")
!      OPEN(UNIT=2,file="bdry.dat")
      !OPEN(UNIT=4,file="n1.dat")
      !OPEN(UNIT=5,file="n2.dat")      
      !open(unit=100,file="bound.dat")
      !OPEN(UNIT=11,file="z2.dat")
      !OPEN(UNIT=12,file="deta3b16h5bound3.dat")
      !open(unit=13,file="0bb2.dat")
 !      OPEN(UNIT=30,file="c0.dat")
 !      OPEN(UNIT=31,file="c1.dat")  
 !      OPEN(UNIT=32,file="c2.dat")  
 !      OPEN(UNIT=33,file="c3.dat")  
 !      OPEN(UNIT=34,file="c4.dat")
      !OPEN(UNIT=35,file="bsita.dat")
      !OPEN(UNIT=36,file="q.dat")
  !     OPEN(UNIT=37,file="dd.dat")
      !OPEN(UNIT=38,file="f.dat")
      !OPEN(UNIT=13,file="WEI.dat")
      !OPEN(UNIT=14,file="WEIDER.dat")
      !OPEN(UNIT=21,file="066ideal.dat")
   !    OPEN(UNIT=301,file="c01.dat")
    !   OPEN(UNIT=302,file="c02.dat")
     !  OPEN(UNIT=303,file="c03.dat")
      ! OPEN(UNIT=304,file="c11.dat")
       !OPEN(UNIT=305,file="c12.dat")
       !OPEN(UNIT=306,file="c13.dat")
       !OPEN(UNIT=307,file="c21.dat")
       !OPEN(UNIT=308,file="c22.dat")
       !OPEN(UNIT=309,file="c23.dat")
   z2h=(0.0,0.0)             !z(2) at a-h
   z1h=(0.0,0.0)
   wei=(0.0,0.0)
   weider=(0.0,0.0)
   a=1.0_8
   k=-1.0/r0 !1.0*0.6/1.67            !0.87=2*1.0*0.67/1.67                                                                         !shooting original value
   d=0.01_8 !d=0.01
   b=1.1_8 !1.2
   deta=1.0e6 !5e4      10e6          !8
   bt=0.8_8
   !wii=-0.04_8
   !h=0.0005_8
   !v0=0.2_8
   !yita=0.1_8 !01
   !dakai=0.01_8 
    qq=1.6  ! the safety fator at the plasma eadge
    m=2   !the number of the mode 
    dw=1.0e-4_8
    dwr=1.0e-4_8
    dv=0.001 

  bn=0

if(qq==1)then
    da=0.03 !the difference of q safe factor
    else
    da=0.01
  endif                 !dw=1.0e-3
    wii=0.07_8
    dakai1=0.0
    v1=0.0
    k11=0.0
  !yita1=0.0
    v0=bn*0.01 !v
    v0=-0.0
   yita=0.2 !eta  !0.1
    dakai=0.10000  !chi 0.3 !0.01
    aa=(2.0*bt)/(qq*r0) !for bita
 !alph=(0.5,-0.01)!save
 alph=(0.02,-0.2) !seems to be ok to repeat the results in Fig.6
 !aa=0.22
! bb2u=(1.,1.0)
 !wiu=1.0

  if(qq/=1)then
  !jkn=(0.16/(qq-1)-0.16/qq)/da
jkn=1.0/da
else
jkn=1.0/da
  !jkn=(1.16-0.16/(qq))/da
endif

!do while(dakai<0.102)
!w1=0.0
!w2=0.0
v0=0.0

do while(v0>=-0.07)
check=0
beta=(0.0,0.0)
ws=(0.0,0.0) !for the interation
w=cmplx(-k*v0+0.00,0.00)
do jk=1,nn
     rr=0.01
     h=hmax
    if(m<=2)then
     z(1)=cmplx(0.1*rr**2,0.1*rr**2) !0.1  !origin value 
     z(2)=cmplx(0.2*rr**1,0.2*rr**1) !0.2  !
    else 
     z(1)=cmplx(rr**(m-1),0.0)
     z(2)=cmplx(m*rr**(m-2),0.0)
    endif  
!       if(abs(real(bb2))>1.0e-8)then
!	apha=(aimag(bb2)-aimag(bb1))/dw
!	dw=-aimag(bb2)/apha
!	bb1=bb2
        ! -k*v0=-k*v0+0.001
!       wii=wii-dw
  !      w=cmplx(-k*v0+wr,wii)
        call rk2(m,aa,w,check,yita,bsita,bt,h,k,Z,k1,k2,k3,k4,k5,k6,rr,z2h,z22h,z1h,dakai)
        call weii(m,aa,w,wr,dakai,yita,bt,h,k,z2h,z22h,z1h,z(1),z(2),wei,weider,wei0,wei00)
        call bbb(z(1),aa,z(2),z2h,z1h,m,w,wii,wr,k,deta,h,a,b,bt,d,bb2,wei,wei0,wei00,weider,yita,dakai,v0)
  
      write(3,*) jk,v0,aimag(bb2),real(bb2)
     if (abs(real(bb2))<er.and.abs(aimag(bb2))<er)then
      exit
      else
      wv(jk)=w
      w=w+alph*bb2
 
    end if
     ! write(2,*) v0,real(w)
    
enddo
   

      
     if (abs(real(bb2))>er.and.abs(aimag(bb2))>er)then
do jj=1,nn-2
      beta=(wv(jj+2)+wv(jj+1))/(wv(jj+1)-wv(jj))+beta
      ws=ws+(wv(jj+2)*wv(jj)-wv(jj+1)**2)/(wv(jj+2)+wv(jj)-2.0*(wv(jj+1)))
enddo
      
     !if (abs(real(bb2))>er.and.abs(aimag(bb2))>er)then
     beta=beta/(nn-2)*1.0
      ws=ws/(nn-2)*1.0
      w=ws
      write(1,*) v0,aimag(w)*r0,(real(w)+k*v0)*r0          
      write(2,*) v0,real(bb2),aimag(bb2)
      
endif
!v1=v0+dv
!w2=aimag(w) ! w2 v0 ;w1 v0-dv
!wou=w1*w2*10e6
!if(wou<0.0)then
!vc=(w2*v1-w1*v0)/(w2-w1)
!write(3,*) dakai,vc

!goto 1
!endif
v0=v0-dv
!w1=aimag(w)
enddo
!enddo
END PROGRAM
 

SUBROUTINE rk2(m,aa,w,ck,yita,bsita,bt,h,k,Z,k1,k2,k3,k4,k5,k6,rr,z2h,z22h,z1h,dakai)
 	USE DATA_MODULE
	IMPLICIT NONE
	INTEGER(4):: i,m,ii
	COMPLEX(8):: w,z2h,z1h,z22h
	real(8):: ck,dakai,yita,bsita,bt,h,k,zz,r,rr,cigema,rtol,f,aa
    COMPLEX(8),DIMENSION(N):: Z,K1,K2,K3,K4,k5,k6
	i=1
	k1=(0.0,0.0)
	k2=(0.0,0.0)
	k3=(0.0,0.0)
	k4=(0.0,0.0)
	ii=((100-rr*100)/(h*100)+1)
        r=rr
	h=hmax
	DO i=1,ii                 !while(r<=1.0) !991 ! 990次 就是991    
	  r=rr
	   ! if(abs(r-0.738)<=0.003)then
	  ! r=0.741
	  ! endif
	  bsita=aa*(1.0-(1.0-r**2)**(miu+1.0))/(r*(miu+1.0)) 
	  f=k*bt+m*bsita/r
      ! write(38,*) r,f
	 
	 CALL RKQC(m,aa,w,dakai,yita,bsita,bt,h,k,z,r,K1,ck)
	 CALL RKQC(m,aa,w,dakai,yita,bsita,bt,h,k,z+0.5_8*K1,r+0.5_8*H,K2,ck)
	 CALL RKQC(m,aa,w,dakai,yita,bsita,bt,h,k,z+0.5_8*K2,r+0.5_8*H,K3,ck)
	 CALL RKQC(m,aa,w,dakai,yita,bsita,bt,h,k,z+1.0_8*K3,r+1.0_8*H,K4,ck)
	 Z=Z+(K1+2.0_8*K2+2.0_8*K3+K4)/6
     
   if(i==ii-2)then
	 z22h=z(2)
	 endif
	 if(i==ii-1)then
	 z2h=z(2)
 	 z1h=z(1)
	endif
!	zz=sqrt(real(z(1))**2+aimag(z(1))**2)
!	 !,aimag(z(1))!最好把在z（1）的模还是实部呢？
	rr=rr+h
	ENDDO 
!	WRITE(10, *) aimag(w),real(z(1)),aimag(z(1))
    END SUBROUTINE rk2

SUBROUTINE RKQC(m,aa,w,dakai,yita,bsita,bt,h,k,Z,r,BACK,check)! 经典的 RungEKuta算法
    USE DATA_MODULE
	IMPLICIT NONE
	INTEGER(4):: m      
	COMPLEX(8):: c0,c1,c2,c3,c4,C00,C11,C22,n1,n2,n3,n4,n5,dd,w,b2,a2
	COMPLEX(8):: c0_1,c0_2,c0_3,c1_1,c1_2,c1_3,c2_1,c2_2,c2_3                  
    complex(8),DIMENSION(N)::Z
	complex(8),DIMENSION(N)::BACK
	real(8):: check,f,ff,g,gg,a1,k02,bsita,r,bt,h,yita,dakai,k,wa,bsitah,bsitahh,bsita2h,bsitader,bsitarder
    real(8):: fh,gh,a1h,k02h,q,pder,fder,frder,aa,h3
	complex(8)::a1rder,a1rderh,a2h,a1r2der,a1der,a12der,a2rder,a2der,a22der,b1der,b12der,b2der,b22der,a23der,b23der,k6,k61,k61f,k612d
     complex(8):: gder,k61k43,jg,a12derh,k42,k62,a1derh,k62d
     
     !  h3=h
      ! h=hmin 

      bsita=aa*(1.0-(1.0-r**2)**(miu+1.0))/(r*(miu+1.0)) !miu current profile
	 !bsita=bsita/bb0
     
	 bsitah=aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/((r-h)*(miu+1.0)) !bsita(r-h)
     !bsitah=bsitah/bb0

	 bsitahh=aa*(1.0-(1.0-(r+h)**2)**(miu+1.0))/((r+h)*(miu+1.0)) !bsita(r+h)
     !bsitahh=bsitahh/bb0

     bsita2h=aa*(1.0-(1.0-(r-2.0*h)**2)**(miu+1.0))/((r-2.0*h)*(miu+1.0))  !bsita(r-2h)
     !bsita2h=bsita2h/bb0
	  
     bsitader=(aa*(1.0-(1.0-r**2)**(miu+1.0))/(r*(miu+1))-aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/((r-h)*(miu+1.0)))/h !derivative bsita
     !bsitader=bsitader/bb0 
	 
	 bsitarder=(aa*(1.0-(1.0-r**2)**(miu+1.0))/((r**2)*(miu+1.0))-aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/(((r-h)**2)*(miu+1.0)))/h !derivative bsita/r
    ! bsitarder=bsitarder/bb0

	 a1=k*r/(r*k*bt+bsita*m)
	
	 f=k*bt+m*bsita/r
	 !f=1-r**2
	
	 !IF(abs(f)<=1e-4)then
	 !r0=r
	 !endif
     pder=((3.0-(1.0+8.0*1**2)**0.5)/2-((3.0-(1.0+8.0*(1-h)**2)**0.5)/2))/h

	 b2=j*bt/(r*f)
    
	 pder=((3.0-(1.0+8.0*r**2)**0.5)/2-((3.0-(1.0+8.0*(r-h)**2)**0.5)/2))/h
	
	 a1h=k*(r-h)/((r-h)*k*bt+bsitah*m) !a1(r-h)
	 
     fh=k*bt+m*bsitah/(r-h) !f(r-h)

     ff=k*bt-m*bsita/r

     g=m*bt/r-k*bsita
	 
	 gh=m*bt/(r-h)-k*bsitah !g(r-h)
     
	 gg=m*bt/r+k*bsita
     
	 a2=j*bsita/(r*f) !complex number
	 
	 a2h=j*bsitah/((r-h)*fh) !a2(r-h)
     
	 k02=k**2+(m/r)**2
	 
	 k02h=k**2+(m/(r-h))**2 ! k02(r-h)
	 
	 wa=f

	 q=r*bt/(bsita*r0)

   gder=(m*bt/r-k*bsita-m*bt/(r-h)+k*bsitah)/h

   fder=(fh-f)/(h*fh) !1/f *F

   frder=(fh*(r-h)-f*r)/(r*fh*(r-h)*h) !derivative 1/fr *f 

   a1rder=(k/(r*k*bt+bsita*m)-k/((r-h)*k*bt+bsitah*m))/h !derivative A1/r

   a1rderh=(k/((r-h)*k*bt+bsitah*m)-k/((r-2.0*h)*k*bt+bsita2h*m))/h !derivative A1/(r-h)

   a1r2der=(k/((r+h)*k*bt+bsitahh*m)+k/((r-h)*k*bt+bsitah*m)-2.0*k/(r*k*bt+bsita*m))/h**2 !2 order derivative A1/r

   a1der=(k*r/(r*k*bt+bsita*m)-k*(r-h)/((r-h)*k*bt+bsitah*m))/h !derivative A1

   a12der=(k*(r+h)/((r+h)*k*bt+bsitahh*m)+k*(r-h)/((r-h)*k*bt+bsitah*m)-2.0*k*r/(r*k*bt+bsita*m))/h**2 !2 order derivative A1

   a1derh=(-k*(r-2.0*h)/((r-2.0*h)*k*bt+bsita2h*m)+k*(r-h)/((r-h)*k*bt+bsitah*m))/h !a1 at r-h

   a12derh=(k*r/(r*k*bt+bsita*m)+k*(r-2.0*h)/((r-2.0*h)*k*bt+bsita2h*m)-2.0*k*(r-h)/((r-h)*k*bt+bsitah*m))/h**2        !2 order derivative a1(r-h)

   a2rder=(j*bsita/(r**2*k*bt+r*m*bsita)-j*bsitah/((r-h)**2*k*bt+(r-h)*m*bsitah))/h !derivative A2/r

  a2der=(j*bsita/(r*k*bt+m*bsita)-j*bsitah/((r-h)*k*bt+m*bsitah))/h  !derivative A2

 a22der=(j*bsitahh/((r+h)*k*bt+m*bsitahh)+j*bsitah/((r-h)*k*bt+m*bsitah)-2.0*j*bsita/(r*k*bt+m*bsita))/h**2 !2 order derivative A2

 b1der=-(m/(k*r*bt+m*bsita)-m/(k*(r-h)*bt+m*bsitah))/h !derivative B1

 b12der=(2.0*m/(k*r*bt+m*bsita)-m/(k*(r-h)*bt+m*bsitah)-m/(k*(r+h)*bt+m*bsitahh))/h**2 !2 order derivative B1

 b2der=(j*bt/(k*r*bt+m*bsita)-j*bt/(k*(r-h)*bt+m*bsitah))/h !!derivative B2

 b22der=(j*bt/(k*(r+h)*bt+m*bsitahh)+j*bt/(k*(r-h)*bt+m*bsitah)-2.0*j*bt/(k*r*bt+m*bsita))/h**2 !2 order derivative B2

 a23der=(j*bsitahh/((r+h)*k*bt+m*bsitahh)+3.0*j*bsitah/((r-h)*k*bt+m*bsitah)-3.0*j*bsita/(r*k*bt+m*bsita)-j*bsita2h/(k*(r-2.0*h)*bt+m*bsita2h))/h**3 !3 order derivative A2

 b23der=(j*bt/(k*(r+h)*bt+m*bsitahh)+3.0*j*bt/(k*(r-h)*bt+m*bsitah)-3.0*j*bt/(k*r*bt+m*bsita)-j*bt/(k*(r-2.0*h)*bt+m*bsita2h))/h**3 !3 order derivative B2

   k61=j*(w*dakai*f*m*2.0*a1*a1rder/(r*k)+j*(g)*(-k02+j*g**2*w*yita-2.0*j*f**2*w*yita)&
      -w*dakai*f*m*2.0*a1h*a1rderh/((r-h)*k)-j*(gh*f/fh)*(-k02h+j*gh**2*w*yita-2.0*j*fh**2*w*yita))/h

   k6=-j*((w*dakai*m*2.0*f)*(a1*a1rder/(r*k))+j*(g)*(-k02+j*g**2*w*yita-2.0*j*f**2*w*yita)) !*2/r**3+a1*a12der/r)
 !K6=K6*-1
   k61f=(2.0*f*w*yita*bsita**2*f/r**2+2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2)+2.0*w*dakai*m**2/(r**4)-&
       (2.0*f*fh*w*yita*bsitah**2/(r-h)**2+2.0*j*gh*bsitah*k*f/(fh*(r-h)**2)-(gh**2*w*yita*bsitah**2*f)/(fh*(r-h)**2)+(2.0*w*dakai*m**2*f)/(fh*(r-h)**4)))/h

   k612d=(j*w**2/(r)-2.0*f**2*w*yita/r+j*g**2/(r)+g**2*w*yita/(r)+j*w*dakai*(j*k02/r+m*a2/r**3)-&
       (j*w**2*f/(fh*(r-h))-2.0*fh*f*w*yita/(r-h)+j*gh**2*f/(fh*(r-h))+(gh**2*w*yita*f)/(fh*(r-h))+j*w*dakai*(j*k02h/(r-h)+m*a2h/(r-h)**3)*f/fh))/h

   k42=j*k6/r-2.0*j*f*((bsita*a1der+bt*b1der)*(w**2-k02+j*g*g*w*yita)+(bt*a1der-bsita*b1der)*(g*w**2/f+2.0*j*g*f*w*yita)&
      -j*w*dakai*(a12der*b1der-a1der*b12der-k02*m*a1**2/(k*r**2)-a1*b1der/r**2))

   k62=-2.0*w*dakai*m*a1/(r**2)        !-2*w*dakai*m*(2*a1der/r-a1/r**2)/f

   k62d=(-(2.0*w*dakai*m*a1)/(r**2)+(2.0*w*dakai*m*a1h*f)/(((r-h)**2)*fh))/h                                  !-2*w*dakai*m*((2*a1der/r-a1/r**2)/f-(2*a1derh/(r-h)-a1h/(r-h)**2)/fh)/h
   
   jg=(-j*g-g*w*yita)*r**4+j*(2.0*w*dakai*m*a1*(2.0*j*k*bsita*r+r**2*gder*(-j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai))/(k6) !*r**4
 !jg=-j*g-g*w*yita+j*(2.0*w*dakai*m*a1*(-2.0*j*k*bsita*r+r**2*gder*(j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai))/(r**4*f*k6)
 !jg=(-j*g-g*w*yita+(2*w*dakai*m*a1*(-2*j*k*bsita*r+r**2*gder*(j-w*yita)-3*g*w*yita*bsita**2*r-2*m*a1*w*dakai))&
 !    /(2*w*dakai*m*f*r*a1**2/(k)+j*k02*g*r**4+g**3*w*yita*r**4-2*f**2*r**4*g*w*yita))
 !chi =-chi 4.21
 
   dd=k42*k62-k6**2+j*k61*(-k62)-j*k6*k62d !*f**2  约掉 (\rho w dakai)**2
   dd=dd*r**6
 
 !dd=((w*dakai*m*2*a1/r)*(a1rder/k)+j*(g/f)*(-k02+j*g**2*w*yita-2*j*f**2*w*yita)&
  !   *(2*w*dakai*m)*((a1der/r**2-2*a1/r**3)/f)+j*(g/f)*(-k02+j*g**2*w*yita-2*j*f**2*w*yita))&
!	 +(4*j*w*dakai*m*a1**2/r**2*k)*((g*a1der/k+m*bt*a1/(r**2*k))*(w**2-k02+j*g**2*w*yita)&
!	 +(g*a1der/k-m*bsita*a1*g/(f*r**2*k))*(w**2+2*j*f**2*w*yita)-j*w*dakai*(a12der*b1der-a1der*b12der-k02*m*a1**2/(k*r**2)&
!	 -a1**2*m/(r**4*k)+a1*a1der*m/(r**3*k)))&
!	 +2*(w*dakai*m*2*a1/r)*(a1rder/k)*((w*dakai*m*2*a1/r)*(a1rder/k)+j*(g/f)*(-k02+j*g**2*w*yita-2*j*f**2*w*yita))
	 !
    k62=-2.0*m*a1/(r**2) 

 
   ! h=h3



 n5=k6*w*dakai/(r*f)

 n4=k61*w*dakai/(r*f)+k6*(-j*w*dakai*(j/r**2+2.0*m*a2der/r+2.0*k*b2der)/f+w*dakai*frder)&
   +k62*((-w**3*dakai-&
    j*k02*w**2*dakai**2+w*dakai*k02-j*g**2*w**2*dakai*yita)/(f*r)+w**2*dakai**2*b2*(a1rder+a12der)&
	-w**2*dakai**2*a2*(b1der/r+b12der))

 n3=-k61*j*w*dakai*(j/r**2+2.0*m*a2der/r+2.0*k*b2der)+k6*((j*w**2/(r)-2.0*f*f*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))&
   +j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2/r+2.0*a2der)*fder/(r)-j*w*dakai*m*((-2.0*a2)/r**2+3.0*a22der)/(r)&
   -j*w*dakai*k*(b2/r+2.0*b2der)*fder-j*w*dakai*k*((2.0*r*b2der-b2)/r**2+3.0*b22der))&
   +k62*((j*w**3*dakai-k02*w**2*dakai**2)&
   *(m*a2/r+2.0*m*a2der+k*(b2+2.0*b2der*r))/(r)-f*(j*w*dakai*k02+g**2*w**2*dakai*yita)*(a2*bsita/r+2.0*bsita*a2der+bt*(b2/r+2.0*b2der))&
   -2.0*f*f*g*w**2*yita*dakai*(2.0*bt*a2der-2.0*bsita*b2der)+(w*dakai)**2*f*(a1rder+a12der)*(b2/r+2.0*b2der)&
   -(w*dakai)**2*f*(b1der/r+b12der)*(a2/r+2.0*a2der))
 
 


 n2=k61*((j*w**2/(r)-2.0*f*f*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))+j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2der/r+a22der)/(r)&
   -j*w*dakai*k*(b2der/r+b22der))&
   +k6*(2.0*f**2*w*yita*bsita**2/r**2+2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2)+2.0*w*dakai*m**2/(r**4)+k612d-&
   (j*w*dakai*m*frder)*(a2der/r+a22der)-(j*w*dakai*m/(r))*((r*a22der-a2der)/r**2+a23der)-(j*w*dakai*k*fder)*& !fder *f
   (b2der/r+b22der)-(j*w*dakai*k)*((r*b22der-b2der)/r**2+b23der))&
   +k62*((j*w**3*dakai-k02*w**2*dakai**2)*((m/r)*(a2der/r+a22der)+k*(b2der/r+b22der)-j*k02/r-m*a2/(r**3))&
   -((k02*j*w*dakai)/(r))*(j*w**2-2.0*f**2*w*yita+j*g**2+g**2*w*yita)-(j*w*dakai*k02+g**2*w**2*dakai*yita)*(bsita*f*(a2der/r+a22der)&
   +bt*f*(b2der/r+b22der)-j*k02/(r)-bsita*a2*f/(r**2))&
   -2.0*f*f*g*w**2*dakai*yita*(bt*(a2der/r+a22der)-bsita*(b2der/r+b22der)-bt*a2/(r**2))+(-w**2+k02-j*(g**2)*w*yita)&
   *(j*w**2/(r)-2.0*f*f*w*yita/r)+(j*w*dakai/r)*(a1rder+a12der)*(j*bt*w**2-2.0*f*f*bt*w*yita-j*g*f*bsita-g*f*bsita*w*yita)-&
   (2.0*j*f*f*g*w*yita+g*w**2)*(j*g/r+g*w*yita/r)+(w*dakai)**2*(a1rder+a12der)*(b2der*f/r+b22der*f-j*k02*bt/(r))+&
   (j*w*dakai)*(b1der/r+b12der)*(-j*g*f*bt-g*f*bt*w*yita-j*bsita*w**2+2.0*f**2*bsita*w*yita)*(1/r)-&
	(w*dakai)**2*(b1der/r+b12der)*(a2rder*f+a22der*f-j*k02*bsita/(r))) 
    


                
 
 n1=k61*(2.0*f*f*w*yita*bsita**2/r**2+2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2)+2.0*w*dakai*m**2/(r**4))&
   +k6*k61f+k62*(-2.0*w**3*dakai*m**2/(r**4)-j*w*dakai*k02*(2.0*f*f*yita*bsita**2/r**2+&
   2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2))-(2.0*j*f*f*g*w*yita+g*w**2)*(2.0*j*bsita*k/r**2-g*w*yita*bsita**2/r**2)&
   -2.0*j*k02*(w*dakai*m)**2/(r**4)-(j*w*dakai*k02+g**2*w**2*dakai*yita)*2.0*j*m*bsita*f/r**3-4.0*j*f**2*w**2*dakai*yita*g*m*bt/(r**3)&
   +(j*w*dakai*f/r**2)*(a1rder+a12der)*(2.0*f*w*yita*(bsita**2)*bt-2.0*j*bsita**2*k+g*w*yita*bsita**3)&
   +(-w**2+k02-j*(g**2)*w*yita)*(2.0*f*f*w*yita*bsita**2/(r**2))-(j*w*dakai*f/(r**2))*(b1der/r+b12der)*(2.0*j*k*bsita*bt-g*w*yita*bsita**2*bt+2.0*f*w*yita*bsita**3)&
   -2.0*j*w**2*dakai**2*m*f*(b1der+b12der)/(r**3))
	
 dakai=dakai     !--chi--------4 -----6 

 n1=n1*r**2

 n2=n2*r**2

 n3=n3*r**2

 n4=n4*r**2

 n5=n5*r**2
                                                                                                                 !jg=0.0
	
 c4=(n5/dd)*jg

 c3=(n4/dd)*jg&
   -(-2.0*j*k*bsita*r+gder*j*r**2+gder*w*yita*r**2+3.0*g*w*yita*bsita**2*r-2.0*m*a1*w*dakai)*&
	 w*dakai/(r**3*(-k6))

 
  c2=(n3/dd)*jg&
     -(-2.0*j*k*bsita*r+gder*j*r**2+gder*w*yita*r**2+3.0*g*w*yita*bsita**2*r-2.0*m*a1*w*dakai)*&
     (-j)*w*dakai*(j/r**2+2.0*m*a2der/r+2.0*k*b2der)/(r**2*(-k6))& !有两个负号相互抵消了的
     -(1.0-j*w*yita-j*w*dakai)/r
!for check
  c2_1=(n3/dd)*jg
  c2_2=-(-2.0*j*k*bsita*r+gder*j*r**2+gder*w*yita*r**2+3.0*g*w*yita*bsita**2*r-2.0*m*a1*w*dakai)*&
     (-j)*w*dakai*(j/r**2+2.0*m*a2der/r+2.0*k*b2der)/(r**2*(-k6))  
  c2_3=-(1.0-j*w*yita-j*w*dakai)/r
!for check
	
  c1=(n2/dd)*jg&
	 +(2.0*j*k*bsita*r+r**2*gder*(-j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai)&
	 *((j*w**2/(r)-2.0*f*f*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))+j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2der/r+a22der)/(r)&
     -j*w*dakai*k*(b2der/r+b22der))/(r**2*(-k6))&
	 +((1.0-j*w*yita-j*w*yita*bsita**2+3.0*j*w*yita*bsita**2-j*w*dakai*(1.0+2.0*j*m*a2))/r**2)
 
!for check
   c1_1=(n2/dd)*jg
   c1_2=+(2.0*j*k*bsita*r+r**2*gder*(-j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai)&
         *((j*w**2/(r)-2.0*f*f*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))+j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2der/r+a22der)/(r)&
     -j*w*dakai*k*(b2der/r+b22der))/(r**2*(-k6))
   c1_3=((1.0-j*w*yita-j*w*yita*bsita**2+3.0*j*w*yita*bsita**2-j*w*dakai*(1.0+2.0*j*m*a2))/r**2)
!for check

 c0=(n1/dd)*jg&
	 +(2.0*j*k*bsita*r+r**2*gder*(-j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai)&
	 *(2.0*f**2*w*yita*bsita**2*(r**2)+(2.0*j*g*bsita*k-g**2*w*yita*bsita**2)*r**2/(1.0)&
	 +2.0*w*dakai*m**2)/(r**6*(-k6))&
	 +(wa**2-w**2+2.0*bsita**2/r**2-2.0*bsita*bsitader/r+4*bsita*bsitarder&
	  -(2.0*j*w*yita)*bsitarder*bsita-3.0*j*w*yita*bsita**4/(r**2)-j*w*dakai*k02)/r !

!for check
 c0_1=(n1/dd)*jg
 c0_2=(2.0*j*k*bsita*r+r**2*gder*(-j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai)&
         *(2.0*f**2*w*yita*bsita**2*(r**2)+(2.0*j*g*bsita*k-g**2*w*yita*bsita**2)*r**2/(1.0)&
         +2.0*w*dakai*m**2)/(r**6*(-k6))
 c0_3=(wa**2-w**2+2.0*bsita**2/r**2-2.0*bsita*bsitader/r+4*bsita*bsitarder&
          -(2.0*j*w*yita)*bsitarder*bsita-3.0*j*w*yita*bsita**4/(r**2)-j*w*dakai*k02)/r
!for check

c0=c0*r**2

c1=c1*r**2

c2=c2*r**2

c3=c3*r**2

c4=c4*r**2
dakai=dakai


 !write(*,*) check
!k62=-2.0*w*dakai*m*a1/(r**2*f)
IF (check > 1000) THEN
 write(*,*)  r
 write(30,*) r,real(c0),aimag(c0) ! ,c1,c2,c3,c4
 write(31,*) r,real(c1),aimag(c1) !c1 
 write(32,*) r,real(c2),aimag(c2) !c2 
 write(33,*) r,real(c3),aimag(c3) !c1 
 write(34,*) r,real(c4),aimag(c4)
 write(37,*) r,real(dd),aimag(dd)
 write(301,*) r,real(c0_1),aimag(c0_1) 
 write(302,*) r,real(c0_2),aimag(c0_2) 
 write(303,*) r,real(c0_3),aimag(c0_3) 
 write(304,*) r,real(c1_1),aimag(c1_1) 
 write(305,*) r,real(c1_2),aimag(c1_2) 
 write(306,*) r,real(c1_3),aimag(c1_3) 
 write(307,*) r,real(c2_1),aimag(c2_1) 
 write(308,*) r,real(c2_2),aimag(c2_2) 
 write(309,*) r,real(c2_3),aimag(c2_3) 
ENDIF
!write(1,*) r,real(a1) !aimag(a1)
!write(2,*) r,real(jg),aimag(jg)
!write(3,*) r,real(k6),aimag(k6)
!write(4,*) r,real(n1),aimag(n1)
!write(5,*) r,real(n2),aimag(n2)
!write(35,*) r,real(bsita) !,aimag(bsita)
!write(36,*) r,q


BACK(1)=H*Z(2)

BACK(2)=-H*(c0*z(1)+c1*z(2))/c2  !+c2*z(3)+c3*z(4))/c4 

!BACK(3)=H*z(4)
!BACK(4)=-H*(c0*z(1)+c1*z(2)+c2*z(3)+c3*z(4))/c4 
return

END SUBROUTINE RKQC

subroutine weii(m,aa,w,wr,dakai,yita,bt,h,k,z2h,z22h,z1h,z1,z2,wei,weider,wei0,wei00) !how to solve z2der
    USE DATA_MODULE
	IMPLICIT NONE
	INTEGER(4):: m
	COMPLEX(8):: c0,c1,c2,c3,c4,C00,C11,C22,n1,n2,n3,n4,n5,dd,w,b2,a2,wei0,weider0
    !complex(8),DIMENSION(N)::Z
        complex(8),DIMENSION(N)::BACK
	real(8):: f,ff,g,gg,aa,a1,k02,bsita,r,bt,h,yita,dakai,k,wa,bsitah,bsitahh,bsita2h,bsitader,bsitarder
        real(8):: fh,gh,a1h,k02h,q,dakai1,yita1,h3,wr
        complex(8)::a1rder,a1rderh,a2h,a1r2der,a1der,a1derh,a12der,a12derh,a2rder,a2der,a22der,b1der,b12der,b2der,b22der,a23der,b23der,k6,k61,k61f,k612d
        complex(8):: gder,fder,frder,k61k43,jg,wei,weider,z2der,z2,z22,z2h,z22h,z1h,z1,k42,k62,k62d,wei00
     r=1.0
  !h3=h
 !h=hmin
      bsita=aa*(1.0-(1.0-r**2)**(miu+1.0))/(r*(miu+1.0)) !miu current profile
     !bsita=bsita/bb0
     
      bsitah=aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/((r-h)*(miu+1.0)) !bsita(r-h)
     !bsitah=bsitah/bb0

      bsitahh=aa*(1.0-(1.0-(r+h)**2)**(miu+1.0))/((r+h)*(miu+1.0)) !bsita(r+h)
     !bsitahh=bsitahh/bb0

     bsita2h=aa*(1.0-(1.0-(r-2.0*h)**2)**(miu+1.0))/((r-2.0*h)*(miu+1.0))  !bsita(r-2h)
     !bsita2h=bsita2h/bb0
	  
     bsitader=(aa*(1.0-(1.0-r**2)**(miu+1.0))/(r*(miu+1))-aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/((r-h)*(miu+1.0)))/h !derivative bsita
     !bsitader=bsitader/bb0 
	 
	 bsitarder=(aa*(1.0-(1.0-r**2)**(miu+1.0))/((r**2)*(miu+1.0))-aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/(((r-h)**2)*(miu+1.0)))/h !derivative bsita/r
    ! bsitarder=bsitarder/bb0

	 a1=k*r/(r*k*bt+bsita*m)
	
	 f=k*bt+m*bsita/r
	 !f=1-r**2
	
	 !IF(abs(f)<=1e-4)then
	 !r0=r
	 !endif
!     pder=((3.0-(1.0+8.0*1**2)**0.5)/2-((3.0-(1.0+8.0*(1-h)**2)**0.5)/2))/h

	 b2=j*bt/(r*f)
    
	! pder=((3.0-(1.0+8.0*r**2)**0.5)/2-((3.0-(1.0+8.0*(r-h)**2)**0.5)/2))/h
	
	 a1h=k*(r-h)/((r-h)*k*bt+bsitah*m) !a1(r-h)
	 
     fh=k*bt+m*bsitah/(r-h) !f(r-h)

     ff=k*bt-m*bsita/r

     g=m*bt/r-k*bsita
	 
	 gh=m*bt/(r-h)-k*bsitah !g(r-h)
     
	 gg=m*bt/r+k*bsita
     
	 a2=j*bsita/(r*f) !complex number
	 
	 a2h=j*bsitah/((r-h)*fh) !a2(r-h)
     
	 k02=k**2+(m/r)**2
	 
	 k02h=k**2+(m/(r-h))**2 ! k02(r-h)
	 
	 wa=f

	 q=r*bt/(bsita*R0)

   gder=(m*bt/r-k*bsita-m*bt/(r-h)+k*bsitah)/h

   fder=(fh-f)/(h*fh) !1/f *F

   frder=(fh*(r-h)-f*r)/(r*fh*(r-h)*h) !derivative 1/fr *f 

   a1rder=(k/(r*k*bt+bsita*m)-k/((r-h)*k*bt+bsitah*m))/h !derivative A1/r

   a1rderh=(k/((r-h)*k*bt+bsitah*m)-k/((r-2.0*h)*k*bt+bsita2h*m))/h !derivative A1/(r-h)

   a1r2der=(k/((r+h)*k*bt+bsitahh*m)+k/((r-h)*k*bt+bsitah*m)-2.0*k/(r*k*bt+bsita*m))/h**2 !2 order derivative A1/r

   a1der=(k*r/(r*k*bt+bsita*m)-k*(r-h)/((r-h)*k*bt+bsitah*m))/h !derivative A1

   a12der=(k*(r+h)/((r+h)*k*bt+bsitahh*m)+k*(r-h)/((r-h)*k*bt+bsitah*m)-2.0*k*r/(r*k*bt+bsita*m))/h**2 !2 order derivative A1

   a1derh=(-k*(r-2.0*h)/((r-2.0*h)*k*bt+bsita2h*m)+k*(r-h)/((r-h)*k*bt+bsitah*m))/h !a1 at r-h

   a12derh=(k*r/(r*k*bt+bsita*m)+k*(r-2.0*h)/((r-2.0*h)*k*bt+bsita2h*m)-2.0*k*(r-h)/((r-h)*k*bt+bsitah*m))/h**2        !2 order derivative a1(r-h)

   a2rder=(j*bsita/(r**2*k*bt+r*m*bsita)-j*bsitah/((r-h)**2*k*bt+(r-h)*m*bsitah))/h !derivative A2/r

  a2der=(j*bsita/(r*k*bt+m*bsita)-j*bsitah/((r-h)*k*bt+m*bsitah))/h  !derivative A2

 a22der=(j*bsitahh/((r+h)*k*bt+m*bsitahh)+j*bsitah/((r-h)*k*bt+m*bsitah)-2.0*j*bsita/(r*k*bt+m*bsita))/h**2 !2 order derivative A2

 b1der=-(m/(k*r*bt+m*bsita)-m/(k*(r-h)*bt+m*bsitah))/h !derivative B1

 b12der=(2.0*m/(k*r*bt+m*bsita)-m/(k*(r-h)*bt+m*bsitah)-m/(k*(r+h)*bt+m*bsitahh))/h**2 !2 order derivative B1

 b2der=(j*bt/(k*r*bt+m*bsita)-j*bt/(k*(r-h)*bt+m*bsitah))/h !!derivative B2

 b22der=(j*bt/(k*(r+h)*bt+m*bsitahh)+j*bt/(k*(r-h)*bt+m*bsitah)-2.0*j*bt/(k*r*bt+m*bsita))/h**2 !2 order derivative B2

  a23der=(j*bsitahh/((r+h)*k*bt+m*bsitahh)+3.0*j*bsitah/((r-h)*k*bt+m*bsitah)-3.0*j*bsita/(r*k*bt+m*bsita)-j*bsita2h/(k*(r-2.0*h)*bt+m*bsita2h))/h**3 !3 order derivative A2

  b23der=(j*bt/(k*(r+h)*bt+m*bsitahh)+3.0*j*bt/(k*(r-h)*bt+m*bsitah)-3.0*j*bt/(k*r*bt+m*bsita)-j*bt/(k*(r-2.0*h)*bt+m*bsita2h))/h**3 !3 order derivative B2

  k61=j*(w*dakai*f*m*2.0*a1*a1rder/(r*k)+j*(g)*(-k02+j*g**2*w*yita-2.0*j*f**2*w*yita)&
      -w*dakai*f*m*2.0*a1h*a1rderh/((r-h)*k)-j*(gh*f/fh)*(-k02h+j*gh**2*w*yita-2.0*j*fh**2*w*yita))/h

   k6=-j*((w*dakai*m*2.0*f)*(a1*a1rder/(r*k))+j*(g)*(-k02+j*g**2*w*yita-2.0*j*f**2*w*yita)) !*2/r**3+a1*a12der/r)
 !K6=K6*-1
   k61f=(2.0*f*w*yita*bsita**2*f/r**2+2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2)+2.0*w*dakai*m**2/(r**4)-&
       (2.0*f*fh*w*yita*bsitah**2/(r-h)**2+2.0*j*gh*bsitah*k*f/(fh*(r-h)**2)-(gh**2*w*yita*bsitah**2*f)/(fh*(r-h)**2)+(2.0*w*dakai*m**2*f)/(fh*(r-h)**4)))/h

   k612d=(j*w**2/(r)-2.0*f**2*w*yita/r+j*g**2/(r)+g**2*w*yita/(r)+j*w*dakai*(j*k02/r+m*a2/r**3)-&
       (j*w**2*f/(fh*(r-h))-2.0*fh*f*w*yita/(r-h)+j*gh**2*f/(fh*(r-h))+(gh**2*w*yita*f)/(fh*(r-h))+j*w*dakai*(j*k02h/(r-h)+m*a2h/(r-h)**3)*f/fh))/h

   k42=j*k6/r-2.0*j*f*((bsita*a1der+bt*b1der)*(w**2-k02+j*g*g*w*yita)+(bt*a1der-bsita*b1der)*(g*w**2/f+2.0*j*g*f*w*yita)&
      -j*w*dakai*(a12der*b1der-a1der*b12der-k02*m*a1**2/(k*r**2)-a1*b1der/r**2))

   k62=-2.0*w*dakai*m*a1/(r**2)        !-2*w*dakai*m*(2*a1der/r-a1/r**2)/f

   k62d=(-2.0*w*dakai*m*a1/(r**2)+2.0*w*dakai*m*a1h*f/(((r-h)**2)*fh))/h                                  !-2*w*dakai*m*((2*a1der/r-a1/r**2)/f-(2*a1derh/(r-h)-a1h/(r-h)**2)/fh)/h
   
   jg=(-j*g-g*w*yita)*r**4+j*(2.0*w*dakai*m*a1*(2.0*j*k*bsita*r+r**2*gder*(-j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai))/(k6)
 !jg=-j*g-g*w*yita+j*(2.0*w*dakai*m*a1*(-2.0*j*k*bsita*r+r**2*gder*(j-w*yita)-3.0*g*w*yita*bsita**2*r+2.0*m*a1*w*dakai))/(r**4*f*k6)
 !jg=(-j*g-g*w*yita+(2*w*dakai*m*a1*(-2*j*k*bsita*r+r**2*gder*(j-w*yita)-3*g*w*yita*bsita**2*r-2*m*a1*w*dakai))&
 !    /(2*w*dakai*m*f*r*a1**2/(k)+j*k02*g*r**4+g**3*w*yita*r**4-2*f**2*r**4*g*w*yita))
 !chi =-chi 4.21
 
   dd=k42*k62-k6**2+j*k61*(-k62)-j*k6*k62d !约掉 (\rho w dakai)**2
   dd=dd*r**6
 
 !dd=((w*dakai*m*2*a1/r)*(a1rder/k)+j*(g/f)*(-k02+j*g**2*w*yita-2*j*f**2*w*yita)&
  !   *(2*w*dakai*m)*((a1der/r**2-2*a1/r**3)/f)+j*(g/f)*(-k02+j*g**2*w*yita-2*j*f**2*w*yita))&
!	 +(4*j*w*dakai*m*a1**2/r**2*k)*((g*a1der/k+m*bt*a1/(r**2*k))*(w**2-k02+j*g**2*w*yita)&
!	 +(g*a1der/k-m*bsita*a1*g/(f*r**2*k))*(w**2+2*j*f**2*w*yita)-j*w*dakai*(a12der*b1der-a1der*b12der-k02*m*a1**2/(k*r**2)&
!	 -a1**2*m/(r**4*k)+a1*a1der*m/(r**3*k)))&
!	 +2*(w*dakai*m*2*a1/r)*(a1rder/k)*((w*dakai*m*2*a1/r)*(a1rder/k)+j*(g/f)*(-k02+j*g**2*w*yita-2*j*f**2*w*yita))
	 !
 k62=-2.0*m*a1/(r**2) 


!h=h3



 n5=k6*w*dakai/(r*f)

 n4=k61*w*dakai/(r*f)+k6*(-j*w*dakai*(j/r**2+2.0*m*a2der/r+2.0*k*b2der)/f+w*dakai*frder)&
   +k62*((-w**3*dakai-&
    j*k02*w**2*dakai**2+w*dakai*k02-j*g**2*w**2*dakai*yita)/(f*r)+w**2*dakai**2*b2*(a1rder+a12der)&
	-w**2*dakai**2*a2*(b1der/r+b12der))

 n3=-k61*j*w*dakai*(j/r**2+2.0*m*a2der/r+2.0*k*b2der)+k6*((j*w**2/(r)-2.0*f*f*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))&
   +j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2/r+2.0*a2der)*fder/(r)-j*w*dakai*m*((-2.0*a2)/r**2+3.0*a22der)/(r)&
   -j*w*dakai*k*(b2/r+2.0*b2der)*fder-j*w*dakai*k*((2.0*r*b2der-b2)/r**2+3.0*b22der))&
   +k62*((j*w**3*dakai-k02*w**2*dakai**2)&
   *(m*a2/r+2.0*m*a2der+k*(b2+2.0*b2der*r))/(r)-f*(j*w*dakai*k02+g**2*w**2*dakai*yita)*(a2*bsita/r+2.0*bsita*a2der+bt*(b2/r+2.0*b2der))&
   -2.0*f*f*g*w**2*yita*dakai*(2.0*bt*a2der-2.0*bsita*b2der)+(w*dakai)**2*f*(a1rder+a12der)*(b2/r+2.0*b2der)&
   -(w*dakai)**2*f*(b1der/r+b12der)*(a2/r+2.0*a2der))
 
 


 n2=k61*((j*w**2/(r)-2.0*f*f*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))+j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2der/r+a22der)/(r)&
   -j*w*dakai*k*(b2der/r+b22der))&
   +k6*(2.0*f**2*w*yita*bsita**2/r**2+2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2)+2.0*w*dakai*m**2/(r**4)+k612d-&
   (j*w*dakai*m*frder)*(a2der/r+a22der)-(j*w*dakai*m/(r))*((r*a22der-a2der)/r**2+a23der)-(j*w*dakai*k*fder)*& !fder *f
   (b2der/r+b22der)-(j*w*dakai*k)*((r*b22der-b2der)/r**2+b23der))&
   +k62*((j*w**3*dakai-k02*w**2*dakai**2)*((m/r)*(a2der/r+a22der)+k*(b2der/r+b22der)-j*k02/r-m*a2/(r**3))&
   -((k02*j*w*dakai)/(r))*(j*w**2-2.0*f**2*w*yita+j*g**2+g**2*w*yita)-(j*w*dakai*k02+g**2*w**2*dakai*yita)*(bsita*f*(a2der/r+a22der)&
   +bt*f*(b2der/r+b22der)-j*k02/(r)-bsita*a2*f/(r**2))&
   -2.0*f*f*g*w**2*dakai*yita*(bt*(a2der/r+a22der)-bsita*(b2der/r+b22der)-bt*a2/(r**2))+(-w**2+k02-j*(g**2)*w*yita)&
   *(j*w**2/(r)-2.0*f*f*w*yita/r)+(j*w*dakai/r)*(a1rder+a12der)*(j*bt*w**2-2.0*f*f*bt*w*yita-j*g*f*bsita-g*f*bsita*w*yita)-&
   (2.0*j*f*f*g*w*yita+g*w**2)*(j*g/r+g*w*yita/r)+(w*dakai)**2*(a1rder+a12der)*(b2der*f/r+b22der*f-j*k02*bt/(r))+&
   (j*w*dakai)*(b1der/r+b12der)*(-j*g*f*bt-g*f*bt*w*yita-j*bsita*w**2+2.0*f**2*bsita*w*yita)*(1/r)-&
	(w*dakai)**2*(b1der/r+b12der)*(a2rder*f+a22der*f-j*k02*bsita/(r))) 
    


                
 
 n1=k61*(2.0*f*f*w*yita*bsita**2/r**2+2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2)+2.0*w*dakai*m**2/(r**4))&
   +k6*k61f+k62*(-2.0*w**3*dakai*m**2/(r**4)-j*w*dakai*k02*(2.0*f*f*yita*bsita**2/r**2+&
   2.0*j*g*bsita*k/(r**2)-g**2*w*yita*bsita**2/(r**2))-(2.0*j*f*f*g*w*yita+g*w**2)*(2.0*j*bsita*k/r**2-g*w*yita*bsita**2/r**2)&
   -2.0*j*k02*(w*dakai*m)**2/(r**4)-(j*w*dakai*k02+g**2*w**2*dakai*yita)*2.0*j*m*bsita*f/r**3-4.0*j*f**2*w**2*dakai*yita*g*m*bt/(r**3)&
   +(j*w*dakai*f/r**2)*(a1rder+a12der)*(2.0*f*w*yita*(bsita**2)*bt-2.0*j*bsita**2*k+g*w*yita*bsita**3)&
   +(-w**2+k02-j*(g**2)*w*yita)*(2.0*f*f*w*yita*bsita**2/(r**2))-(j*w*dakai*f/(r**2))*(b1der/r+b12der)*(2.0*j*k*bsita*bt-g*w*yita*bsita**2*bt+2.0*f*w*yita*bsita**3)&
   -2.0*j*w**2*dakai**2*m*f*(b1der+b12der)/(r**3))
	
 dakai=dakai     !--chi--------4 -----6 


 n1=n1*r**6
 n2=n2*r**6
 n3=n3*r**6
 n4=n4*r**6
 n5=n5*r**6
	

 k62=-2.0*w*dakai*m*a1/(r**2)

 z2der=(z2-z2h)/h                                                !derivative a(2) at a
!z2der=0.0
 weider=(n1*z1+n2*z2+n3*z2der)/dd !derivative of  displace in \eta direction  !+n4*(((z2-z2h)/h-(z2h-z22h)/h)/h)
 wei=(2.0*f*f*w*yita*bsita**2*(r**2)+(2.0*j*g*bsita*k-g**2*w*yita*bsita**2)*r**2&
      +2.0*w*dakai*m**2)*z1/(-k6*r**4)+&
      ((j*(w)**2/(r)-2.0*f*f*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))+j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2der/r+a22der)/(r)&
     -j*w*dakai*k*(b2der/r+b22der))*z2/(-k6)&
     +j*w*dakai*(j/r**2+2*m*a2der/r+2*k*b2der)*z2der/(k6)&
     +j*(2*w*dakai*m*a1)*weider/(r**2*k6) ! -w*dakai*(((z2-z2h)/h-(z2h-z22h)/h)/h)/(f*r*k6)
!wei=j*wei


  
 !dakai1=dakai
  !yita1=yita
z22=z2
!z2=0.0
weider0=(n2*z1+n3*z2)/dd
wei00=((j*w**2/(r)-2*f**2*w*yita/r+j*g**2/(r)+g**2*w*yita/(r))+j*w*dakai*(j*k02/r+m*a2/r**3)-j*w*dakai*m*(a2der/r+a22der)/(r)&
     -j*w*dakai*k*(b2der/r+b22der))*z1/(-k6)&
     +j*w*dakai*(j/r**2+2*m*a2der/r+2*k*b2der)*z2/(k6)&
     +j*(2*w*dakai*m*a1)*weider0/(r**2*k6)
z2=z22
!wei0=((-2*f*w*yita/r+g**2*w*yita/(f*r))+j*w*dakai*(j*k02/r+m*a2/r**3)/f-j*w*dakai*m*(a2der/r+a22der)/(r*f)&
 !    -j*w*dakai*k*(b2der/r+b22der)/f)*z1/(-k6)&
   !  +j*w*dakai*(j/r**2+2*m*a2der/r+2*k*b2der)*z2/(f*k6)&
   !  +j*(2*w*dakai*m*a1)*weider0/(r**2*f*k6)

!



	 
!	 dakai=0.0_8
!	 yita=0.0_8
     




!	 dakai=dakai1
!	 yita=yita1

! wei=(((-2*j*k*bsita*r+r**2*gder*(j-w*yita)-3*g*w*yita*bsita**2*r-2*m*a1*w*dakai)&
! /(r**6*k6))&
! (-2*j*k*bsita*r+r**2*gder*(j-w*yita)-3*g*w*yita*bsita**2*r-2*m*a1*w*dakai)&
! *((j*w**2/(f*r)-2*j*f*w*yita/r+j*g**2/(f*r)+g**2*w*yita/(f*r))+j*w*dakai*(j*k02/r+m*a2/r**3)/f-j*w*dakai*m*(a2der/r+a22der)/(r*f)&
! -j*w*dakai*k*(b2der/r+b22der)/f)*z2/(r**2*k6)&
! -(2*j*k*bsita*r-gder*j*r**2+gder*w*yita*r**2+3*g*w*yita*bsita**2*r+2*m*a1*w*dakai)&
! *j*w*dakai*(j/r**2+2*m*a2der/r+2*k*b2der)*z2der/(f*r**3*k6)&
! -j*(2*w*dakai*m*a1)/(r**2*f*k6)*weider !displace in \eta direction
! WRITE(13, *) aimag(w),real(WEI),aimag(WEI)
! WRITE(14, *) aimag(w),real(WEIDER),aimag(WEIDER)

end subroutine weii



subroutine bbb(z1,aa,z2,z2h,z1h,m,w,wii,wr,k,deta,h,a,b,bt,d,bb,wei,wei0,wei00,weider,yita,dakai,v0)
 USE DATA_MODULE
 IMPLICIT NONE
 REAL(8):: wi,deta,h,h2,aa,a,b,d,kbder,kb,ibder,ib2der,kb2der,ib,kader,ka,iader,ia,kb2,ka2,k,mm34,v0 !miu magnetic conduction
 integer m
 complex(8):: bb,z1,z2,w,w1,lam,mm1,mm,wei,wei0,wei00,weider,z2h,z2der,cc2,c112,b1tp,b1sitap,z1der,z1h,p1,pder,s1,a2,bb2
 external bessi,bessi1,bessk,bessk1
 real(8)::bessi,bessi1,bessk,bessk1,beta,bsitader,f,bsitah,bt,bsita,r,fdeer,ia2der,ka2der,ff,yita,dakai,g,fh,fder,a1,qa,p,wr,wii
 ! sqrt(wi*deta) 不能用变量代替？？
 r=1.0 !r=a

     bsita=aa*(1.0-(1.0-r**2)**(miu+1.0))/(r*(miu+1.0)) !miu current profile
	 !bsita=bsita/bb0
     
	 bsitah=aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/((r-h)*(miu+1.0)) !bsita(r-h)
     !bsitah=bsitah/bb0

	 bsitader=(aa*(1.0-(1.0-r**2)**(miu+1.0))/(r*(miu+1))-aa*(1.0-(1.0-(r-h)**2)**(miu+1.0))/((r-h)*(miu+1.0)))/h 


 qa=bt*r/(bsita*r0)
 a1=k*r/(r*k*bt+bsita*m)
 kb2=m**2/b**2+k**2
 ka2=m**2/a**2+k**2

 f=k*bt+m*bsita/1.0 !f at a
 fh=k*bt+m*bsitah/(r-h)
 !fder=(fh-f)/(h*f*fh)
 
 a2=j*bsita/(r*f)
 g=m*bt/r-k*bsita
 ff=k*bt-m*bsita/1.0 
 fdeer=(m*bsita/1.0-m*bsitah/1.0)/h
 
 
! bsitader=(bsita-bsitah)/h
  h2=h
  h=0.00001
  k=-k
 if(m==1) then
  kader=(bessk1(k*a)-bessk1(k*(a-h)))/(k*h)
  kader=-kader
 
  ka2der=(bessk1(k*a)+bessk1(k*(a-2*h))-2*bessk1(k*(a-h)))/((k*h)**2)

  kbder=(bessk1(k*b)-bessk1(k*(b-h)))/(k*h)
  kbder=-kbder
 
 kb2der=(bessk1(k*b)+bessk1(k*(b-2*h))-2*bessk1(k*(b-h)))/((k*h)**2)
  ka=bessk1(k*a)
  kb=bessk1(k*b)
 else
  kader=(bessk(m,k*a)-bessk(m,k*(a-h)))/(k*h)
  kader=-kader

  ka2der=(bessk(m,k*a)+bessk(m,k*(a-2*h))-2*bessk(m,k*(a-h)))/((k*h)**2)
  kbder=(bessk(m,k*b)-bessk(m,k*(b-h)))/(k*h)
  kbder=-kbder
 
  kb2der=(bessk(m,k*b)+bessk(m,k*(b-2*h))-2*bessk(m,k*(b-h)))/((k*h)**2)
  ka=bessk(m,k*a)
  kb=bessk(m,k*b)
  endif
 
 if(m==1) then
 iader=(bessi1(k*a)-bessi1(k*(a-h)))/(k*h)
 iader=-iader 
 
 ia2der=(bessi1(k*a)+bessi1(k*(a-2*h))-2*bessi1(k*(a-h)))/((k*h)**2)
 ibder=(bessi1(k*b)-bessi1(k*(b-h)))/(k*h)
 ibder=-ibder
 
 ib2der=(bessi1(k*b)+bessi1(k*(b-2*h))-2*bessi1(k*(b-h)))/((k*h)**2)
 ia=bessi1(k*a)
 ib=bessi1(k*b)
 else
 iader=(bessi(m,k*a)-bessi(m,k*(a-h)))/(k*h)
 iader=-iader

 ia2der=(bessi(m,k*a)+bessi(m,k*(a-2*h))-2*bessi(m,k*(a-h)))/((k*h)**2)
 ibder=(bessi(m,k*b)-bessi(m,k*(b-h)))/(k*h)
 ibder=-ibder

 ib2der=(bessi(m,k*b)+bessi(m,k*(b-2*h))-2*bessi(m,k*(b-h)))/((k*h)**2)
 ia=bessi(m,k*a)
 ib=bessi(m,k*b)
 endif
 k=-k 
  h=h2
 !lam=cdsqrt(-j*w*deta)
 !mm=(1-cdexp(2*lam*d))*((k*kbder)**2-(kb2*kb/lam)**2)/((k*kbder-kb2*kb/lam)&
 !*(k*ibder+kb2*ib/lam)*cdexp(2*lam*d)-(k*kbder+kb2*kb/lam)*(k*ibder-kb2*ib/lam))
 !mm1=(cdexp(-2*cdsqrt(-j*w*deta/(0.67**2))*d)-1)*((k*kbder)**2-(kb2*kb/(cdsqrt(-j*w*deta/(0.67**2))))**2)/((k*kbder-kb2*kb/(cdsqrt(-j*w*deta/(0.67**2))))&
 !*(k*ibder+kb2*ib/(cdsqrt(-j*w*deta/(0.67**2))))-cdexp(-2*cdsqrt(-j*w*deta/(0.67**2))*d)*(k*kbder+kb2*kb/(cdsqrt(-j*w*deta/(0.67**2))))*(k*ibder-kb2*ib/(cdsqrt(-j*w*deta/(0.67**2))))) 
 !mm1=(-2*lam*d)*((k*kbder)**2-(kb2*kb/lam)**2)/((k*kbder-kb2*kb/lam)&
 !*(k*ibder+kb2*ib/lam)-(1-2*lam*d)*(k*kbder+kb2*kb/lam)*(k*ibder-kb2*ib/lam))
 !mm1=(cdexp(-2.0*cdsqrt(-j*w*deta)*d)-1.0)*((k*kbder)**2-(kb2*kb/(cdsqrt(-j*w*deta)))**2)/((k*kbder-kb2*kb/(cdsqrt(-j*w*deta)))&
 !*(k*ibder+kb2*ib/(cdsqrt(-j*w*deta)))-cdexp(-2.0*cdsqrt(-j*w*deta)*d)*(k*kbder+kb2*kb/(cdsqrt(-j*w*deta)))*(k*ibder-kb2*ib/(cdsqrt(-j*w*deta)))) 
 !c1/c2
   wr=real(w)
wii=aimag(w)
     w1=w
 ! w=j*aimag(w)
     w=cmplx(wr+k*v0,wii) 
     mm1=-j*w*(kbder**2)/(k*(ibder*kb2der-ib2der*kbder)/(deta*d)+j*w*ibder*kbder)
     w=w1

 !mm1=-j*j*AIMAG(w)*kbder**2/(k*(ibder*kb2der-ib2der*kbder)/(deta*d)+j*j*aimag(w)*ibder*kbder)
 !mm1=1

 z2der=(z2-z2h)/h
 z1der=(z1-z1h)/h
 !z2=z1der
 p=2.0*aa**2*bb0**2/(miu0*(miu+1))
 beta=p*miu0/(bb0**2) !0.5
 !beta=2.0*beta
 !pder=((3.0-(1.0+8.0*1**2)**0.5)/2-((3.0-(1.0+8.0*(1-h)**2)**0.5)/2))/h
  pder=-((r**2-3.0*r**4/4.0+r**6/6.0)-((r-h)**2-0.75*(r-h)**4+(r-h)**6/6.0))/h
  p1=-1.0*beta*(z1/1.0)*pder !z1/1 disturbe displace at r=a beta~3.8

  b1sitap=j*k*wei-z1*bsitader/r-bsita*(z2-z1/r)/r 

  b1tp=-bt*z1/r**2-bt*(z2-z1/r)/r-j*m*wei  !derivative of btp =0

 cc2=j*f*z1/(k*(mm1*iader+kader))  !b contimue in radial direction

 !cc2=(j*fdeer*z1+j*f*(z2-z1))/(ka2*(mm1*ia+ka))!B wu yuan xing bu chengli 
 
 !cc2=0.0
 
 !p1=0.0
  !cc2=(j*f*(z2-z1)+j*fdeer*z1)/(ka2*(mm1*ia+ka))
 
 s1=-g*w*wei+j*w*z2-j*w*bsita**2*z1
 
 !s1=0
!p1=0.0
!c112=j*f*z1/(abs(k)*kader)
!c112=j*f*z1/(k*kader)
!bb=(p1+bsita*(b1sitap)+bt*(b1tp)& !P1 should be complex number when there is rotation !p1 is real number when there is not rotation
 !   -(bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))) !

!bb2=(p1+bsita*(b1sitap)+bt*(b1tp)& !P1 should be complex number when there is rotation !p1 is real number when there is not rotation
 !   -(bsita*j*m*c112*(ka)/r+bt*j*k*c112*(ka))) !

!bb2= -(bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))
 
!bb2=z2-z1*((m*(m-n*qa)**2)*((1.0+(a/b)**(2.0*m))/(1.0-(a/b)**(2.0*m))-2.0/(m-n*qa))/((w**2*a**2/bsita**2)-(m-n*qa)**2))
 !bb=(p1+bsita*(b1sitap)+bt*(b1tp)& !1 should be complex number when there is rotation !p1 is real number when there is not rotation
 !   -((z1/r)*(-(bsita)**2/r-bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))) !+2.0*bsita**2*z1+yita*s1+3.0*yita*bsita**2*j*w*z1+j*dakai*w*(z2)
! p1=0.0
! bb=p1+bsita*(b1sitap)+bt*(b1tp)+yita*s1+3.0*yita*bsita**2*j*w*z1-2.0*bsita**2*z1& !
 !   +j*w*dakai*(z2-z1-2*j*m*a2*z1)-((z1/r)*(-(bsita)**2/r-bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka)) 
 ! old
!p1=0.0
 bb=p1+bsita*(b1sitap)+bt*(b1tp)+yita*s1+3.0*yita*bsita**2*(-g*w*wei00+j*w*z1)& 
  +j*w*dakai*(z2)&
   -((z1/r)*(-(bsita)**2/r-bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))

!nn
 !  bb=p1+bsita*(b1sitap)+bt*(b1tp)+yita*s1+3.0*yita*bsita**2*(-g*w*wei00+j*w*z1)& 
!     +j*w*dakai*(2.0*z2-z1-2.0*j*m*a2*z1-2.0*j*m*a1*wei00+2.0*j*m*(a1*wei+a2*z2))&
  !   -((z1/r)*(-(bsita)**2/r-bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))

!B1
!  bb=p1+bsita*(b1sitap)+bt*(b1tp)+yita*s1+3.0*yita*bsita**2*(-g*w*wei00+j*w*z1)& 
 !    +j*w*dakai*(z2-2.0*j*m*a2*z1+2.0*j*m*(a1*wei+a2*z2))&
  !   -((z1/r)*(-(bsita)**2/r-bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))

!B00  
!bb=p1+bsita*(b1sitap)+bt*(b1tp)+yita*s1+3.0*yita*bsita**2*(-g*w*wei00+j*w*z1)& 
 !    +j*w*dakai*(z2+2.0*j*m*(a1*wei+a2*z2))&
  !   -((z1/r)*(-(bsita)**2/r-bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))


!B2
 !bb=p1+bsita*(b1sitap)+bt*(b1tp)+yita*s1+3.0*yita*bsita**2*(-g*w*wei00+j*w*z1)& 
  !  +j*w*dakai*(z2-2.0*j*m*a2*z1-j*m*a1*wei00+j*m*(a1*wei+a2*z2))&
   !-((z1/r)*(-(bsita)**2/r-bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))


   !  +(-2*j*bsita*k)*wei00-(3*yita*bsita**2*g*w+2*m*w*dakai*a1)*wei0
!  bb=-bb
!
!bb=bsita*(b1sitap)+bt*(b1tp)&                                                                 !
!-((z1/r)*(-(bsita)**2/r+bsita*bsitader)+bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka)) 

 !bb=((-j*w*yita*bsita**2+3.0*j*w*yita*bsita**2-j*w*dakai*(1.0+2.0*j*m*a2)))*z1&
 !+(-1.0+j*w*yita+j*w*dakai)*z2-1.0*j*(fder*r*cc2*k*(mm1*iader+kader)+fder*cc2*k*(mm1*iader+kader)+&
 !fder*cc2*k*(mm1*ia2der+ka2der))+(-j*g-w*yita*g)*wei-(-j*g)*wei0

!bb=p1+bsita*(b1sitap)+bt*(b1tp)-(bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka))
!bb=cc2*k*(mm1*ia2der+ka2der)+j*m*b1sitap+j*k*b1tp !dot B =0
!边界上r=a，都是无量纲化的量
!valuum Bsita bt  which are in the plasma
!bb=a*ka2*(mm1*ia+ka)*z2-z1*(k*a*(mm1*iader+kader)+ka2*(mm1*ia+ka))
!bb=bsita*j*m*cc2*(mm1*ia+ka)/r+bt*j*k*cc2*(mm1*ia+ka)
!bb=-ka2*cc2*(mm1*ia+ka)+k*cc2*(mm1*iader+kader)+k**2*cc2*(mm1*ia2der+ka2der)  !+(j*f*z1+j*fdeer*z1+j*f*(z2-z1))
!bb=-ka2*(mm1*ia+ka)+k*(mm1*iader+kader)+k**2*(mm1*ia2der+ka2der)
!bb=cc2*ka2*(mm1*ia+ka)-(j*f*(z2-z1)+j*fdeer*z1)
!bb=(j*f*z1+j*fdeer*z1+j*f*(z2-z1))/(ka2*(mm1*ia+ka))-((j*f*(z2-z1)+j*fdeer*z1)/(ka2*(mm1*ia+ka)))
!bb=j*f*cc2*(mm1*ia+ka)*ka2+f*ff*z1+f*f*z2
!s1=2*pi*0.67*bsita*bb0*1.0e-6/miu0
! beta=beta*0.67*bb0*miu0*100/(2*pi*0.67*bsita*1.0e-6)
! mm1=((mm1*ia+ka)/((mm1*iader+kader)))
 WRITE(12, *) aimag(w),real(bb),aimag(bb)
 write(13,*)  aimag(w),real(bb2),aimag(bb2)
 end subroutine bbb

FUNCTION bessi1(x) !bessel I 1
REAL(8) bessi1,x
REAL(8) ax
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,&
                 q3,q4,q5,q6,q7,q8,q9,y
SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,&
     0.51498869d0,0.15084934d0,0.2658733d-1,&
	 0.301532d-2,0.32411d-3/
DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,&
     -0.3988024d-1,-0.362018d-2,0.163801d-2,&
	 -0.1031555d-1,0.2282967d-1,-0.2895312d-1,&
     0.1787654d-1,-0.420059d-2/
if (abs(x)<3.75) then
  y=(x/3.75)**2
  bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
  ax=abs(x)
  y=3.75/ax
  bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*&
         (q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
  if(x<0.)bessi1=-bessi1
endif
END FUNCTION bessi1



FUNCTION bessi(n,x) !bessel I N 
INTEGER n,IACC
REAL(8) bessi,x,BIGNO,BIGNI
PARAMETER (IACC=40,BIGNO=1.0e10,BIGNI=1.0e-10)
!USES bessi0
INTEGER j,m
REAL(8) bi,bim,bip,tox,bessi0
if (n<2) pause 'bad argument n in bessi'
if (x==0.) then
  bessi=0.
else
  tox=2.0/abs(x)
  bip=0.0
  bi=1.0
  bessi=0.
  m=2*((n+int(sqrt(float(IACC*n)))))
  do j=m,1,-1
    bim=bip+float(j)*tox*bi
    bip=bi
    bi=bim
    if (abs(bi)>BIGNO) then
      bessi=bessi*BIGNI
      bi=bi*BIGNI
      bip=bip*BIGNI
    endif
    if (j==n) bessi=bip
  end do
  bessi=bessi*bessi0(x)/bi
  if (x<0..and.mod(n,2)==1) bessi=-bessi
endif
END FUNCTION bessi

FUNCTION bessk1(x) !bessel K1
REAL(8) bessk1,x
!USES bessi1
REAL(8) bessi1
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,&
                        q1,q2,q3,q4,q5,q6,q7,y
SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,&
     -0.67278579d0,-0.18156897d0,-0.1919402d-1,&
	 -0.110404d-2,-0.4686d-4/
DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,&
     -0.3655620d-1,0.1504268d-1,-0.780353d-2,&
	 0.325614d-2,-0.68245d-3/
if (x<=2.0) then
  y=x*x/4.0
  bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*&
         (p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
  y=2.0/x
  bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+&
         y*(q5+y*(q6+y*q7))))))
endif
END FUNCTION bessk1

FUNCTION bessk(n,x) !bessel k n
INTEGER n
REAL(8) bessk,x
!USES bessk0,bessk1
INTEGER j
REAL(8) bk,bkm,bkp,tox,bessk0,bessk1
if (n<2) pause 'bad argument n in bessk'
tox=2.0/x
bkm=bessk0(x)
bk=bessk1(x)
do j=1,n-1
  bkp=bkm+j*tox*bk
  bkm=bk
  bk=bkp
end do
bessk=bk
END FUNCTION bessk

FUNCTION bessk0(x)
REAL(8) bessk0,x
!USES bessi0
REAL(8) bessi0
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,&
                        q1,q2,q3,q4,q5,q6,q7,y
SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,&
  0.23069756d0,0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,&
  0.2189568d-1,-0.1062446d-1,0.587872d-2,-0.251540d-2,&
  0.53208d-3/
if (x<=2.0) then
  y=x*x/4.0
  bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+y*&
         (p4+y*(p5+y*(p6+y*p7))))))
else
  y=(2.0/x)
  bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*&
         (q5+y*(q6+y*q7))))))
endif
END FUNCTION bessk0

FUNCTION bessi0(x)
REAL(8) bessi0,x
REAL(8) ax
SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,&
     q7,q8,q9
DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,&
     3.0899424d0,1.2067492d0,0.2659732d0,0.360768d-1,&
	 0.45813d-2/
DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,&
     0.1328592d-1,0.225319d-2,-0.157565d-2,&
	 0.916281d-2,-0.2057706d-1,0.2635537d-1,&
     -0.1647633d-1,0.392377d-2/
if (abs(x)<3.75) then
  y=(x/3.75)**2
  bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
else
  ax=abs(x)
  y=3.75/ax
  bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*&
         (q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
endif
END FUNCTION bessi0

 !FUNCTION bsita
 !FUNCTION b1
! FUNCTION a1
