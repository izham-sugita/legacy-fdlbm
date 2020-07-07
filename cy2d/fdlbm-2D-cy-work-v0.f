!*********************************************************************************
!  2-D Flow Simulation Using Finite Difference Lattice Boltzmann Method          * 
!  Flow past cylinder                                                            *
!				by Muhammad Izham@Sugita                         *  
!	   (code modified from courtesy of Kotera Kishi)                         * 
! Department of Mechanical and System Engineering,                               *
! Kyoto Institute of Technology,Kyoto                                            *                                   
!*********************************************************************************
      implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)   
	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax
        common /A3/TAUINV,ROUINV,ROUCs2,ROUCs2INV
 	common /A4/ dt

        dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)
	dimension udf(iXm,iYm),vdf(iXm,iYm),pdf(iXm,iYm)
	dimension Palph(-1:iXm+2,-1:iYm+2,0:8)
	dimension Ealph(0:iXm+1,0:iYm+1,0:8)

	dimension P1(-1:iXm+2,-1:iYm+2,0:8)
	dimension E1(-1:iXm+2,-1:iYm+2,0:8)        

	dimension w(0:8),ex(0:8),ey(0:8)

!     Contravariant data array
	dimension cv_ksi(iXm,iYm,0:8)
	dimension cv_eta(iXm,iYm,0:8)

!     Grid data array
	dimension x(iXm,iYm)
        dimension y(iXm,iYm)
	dimension yacobi(iXm,iYm)
        dimension x_ksi(iXm,iYm) 
	dimension x_eta(iXm,iYm) 
        dimension y_ksi(iXm,iYm)
	dimension y_eta(iXm,iYm)

      character(len=10) ts,tf

      call READIN(L,x,y)

      call date_and_time(time=ts)
	write(*,*) ts

      call COEFFI(ex,ey,w,cv_ksi, cv_eta, yacobi,
     &  x_ksi, x_eta, y_ksi, y_eta )

	call IC(p,u,v,ex,ey,w,Palph,P1,E1)

      call SOLVE(p,u,v,ex,ey,w,Palph,P1,E1,
     &	udf,vdf,pdf,cv_ksi,cv_eta, x, y, 
     &  yacobi, x_ksi, x_eta, y_ksi, y_eta )
	
      call date_and_time(time=tf)
	write(*,*)ts
      write(*,*)tf

	stop
	end
	
********************************************************
*    Input Data                                        *
********************************************************
     	subroutine READIN(L,x,y)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)   
	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax
	common /A4/ dt

	dimension x(iXm,iYm)
        dimension y(iXm,iYm)

        write(6,*) 'Re=?'
	read(5,*) Re
	write(6,*) 'U0=?'
	read(5,*) U0
	write(*,*) "dt?"
	read(*,*) dt
       	write(*,*) "File number for calculation parameter(3 digit)"
       	write(*,*) "Do not use 100"
	read(*,*) ifile

!	D is density

        D=2.7d0

	write(*,*)"Opening grid data file : JIFCylinder.txt"
	open(11,file="JIFCylinder.txt")
	rewind(11)
	read(11,*) iXmax, iYmax
	do j=1, iYmax
	do i=1, iXmax
	read(11,*) x(i,j), y(i,j)
	end do
	end do
	close(11)

	iXsmax=iXmax-1
	iYsmax=iYmax-1
!	Characteristic length calculation
	diameter = sqrt(x(1,1)**2 + y(1,1)**2)
	diameter = 2.0*diameter

	write(*,*) "Diameter is ", diameter
	write(*,*) "Radius is ", diameter*0.5

	radouter2=x(iXmax-3,1)**2 + y(iXmax-3,1)**2
	diamout=2.0*sqrt(radouter2)
	ratio_r1r2=0.5*diamout/(0.5*diameter)

	write(*,*) "Physical domain radius ratio is ", ratio_r1r2
	write(*,*)"Characteristic length ", diameter
	
	rNYU=diameter*U0/Re

	TAU=3.0d0*rNYU 

	write(*,*) "Tau is ", TAU
	write(*,*) "dt/TAU is ", dt/TAU
        write(*,*)"Residual history is file no. 100"

        write(ifile,*) "Re=",Re
        write(ifile,*) "U0=",U0
        write(ifile,*) "dt=",dt
        write(ifile,*) "tau=",TAU
        write(ifile,*) "Grid :",iXmax, iYmax

	return
	end
********************************************************
*    COEFFICIENT                                       *
********************************************************
      subroutine COEFFI(ex,ey,w,cv_ksi,cv_eta,yacobi,
     &  x_ksi, x_eta, y_ksi, y_eta )
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)   
	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax
	common /A3/TAUINV,ROUINV,ROUCs2,ROUCs2INV


	dimension w(0:8),ex(0:8),ey(0:8)

	dimension yacobi(iXm,iYm)
        dimension x_ksi(iXm,iYm)
	dimension x_eta(iXm,iYm)
        dimension y_ksi(iXm,iYm)
	dimension y_eta(iXm,iYm)

!     Contravariant data array
	dimension cv_ksi(iXm,iYm,0:8)
	dimension cv_eta(iXm,iYm,0:8)

!	dimension cv_ksi(iXmax,iYmax,0:8)
!	dimension cv_eta(iXmax,iYmax,0:8)

	TAUINV=1.0d0/TAU
	write(*,*)"1/Tau=",TAUINV
    
        Cs=1.0d0/dsqrt(3.0d0)
	ROUINV=1.0d0/D
	ROUCs2=D*Cs*Cs
	ROUCs2INV=1.0d0/(D*Cs*Cs)

	open(14,file="JMDCylinder.txt")
	rewind(14)
	read(14,*) itest_x, jtest_y
	if(itest_x .eq. iXmax .and. jtest_y .eq. iYmax) then
	write(*,*)"Opening correct file"
	write(*,*)"imax ",iXmax,"jmax",iYmax
	else
	write(*,*)"Error! Wrong file!"
	stop
	endif
	do j=1, iYmax
	do i=1, iXmax
	read(14,*) yacobi(i,j),x_ksi(i,j),x_eta(i,j),y_ksi(i,j),y_eta(i,j)
	end do
	end do
	close(14)
! 
        ex(0)= 0.0d0
	ex(1)= 1.0d0
	ex(2)= 0.0d0
	ex(3)=-1.0d0
	ex(4)= 0.0d0
	ex(5)= 1.0d0
	ex(6)=-1.0d0
	ex(7)=-1.0d0
	ex(8)= 1.0d0
!
        ey(0)= 0.0d0
	ey(1)= 0.0d0
	ey(2)= 1.0d0
	ey(3)= 0.0d0
	ey(4)=-1.0d0
	ey(5)= 1.0d0
	ey(6)= 1.0d0
	ey(7)=-1.0d0
	ey(8)=-1.0d0
!
        w1=1.0d0/9.0d0
	w2=1.0d0/36.0d0
	w(0)=4.0d0/9.0d0
        w(1)=w1
        w(2)=w1
	w(3)=w1
	w(4)=w1
	w(5)=w2
	w(6)=w2
	w(7)=w2
	w(8)=w2

!	Calculating contravariant velocity
	do k=0,8
	do j=1, iYmax
	do i=1, iXmax
	cv_ksi(i,j,k) = ex(k)*(y_eta(i,j)/yacobi(i,j))
     &	 + ey(k)*(-x_eta(i,j)/yacobi(i,j))

	cv_eta(i,j,k) = ex(k)*(-y_ksi(i,j)/yacobi(i,j))
     &	 + ey(k)*(x_ksi(i,j)/yacobi(i,j))
	end do
	end do
	end do

	return
	end

********************************************************
*    INITIAL CONDITION                                 *
********************************************************
        subroutine IC(p,u,v,ex,ey,w,Palph,P1,E1)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)  
	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax
        common /A3/TAUINV,ROUINV,ROUCs2,ROUCs2INV
	
        dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)
	dimension Palph(-1:iXm+2,-1:iYm+2,0:8)
	dimension P1(-1:iXm+2,-1:iYm+2,0:8)
	dimension E1(0:iXm+1,0:iYm+1,0:8)
	dimension w(0:8),ex(0:8),ey(0:8)

	do j=1,iYmax
        do  i=1,iXmax
        u(i,j)=U0
	v(i,j)=0.0d0
	p(i,j)=ROUCs2
	
	if(i .eq. 1) then
	u(i,j)=0.0
	end if

	end do
	end do

        do 120 k=0,8
	  do 120 j=1,iYmax
	    do 120 i=1,iXmax
	Palph(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

        P1(i,j,k) = Palph(i,j,k)
        E1(i,j,k) = Palph(i,j,k)
	   
  120 continue
      return
	end
********************************************************
*     SOLVE                                            *
********************************************************
       subroutine SOLVE(p,u,v,ex,ey,w,Palph,P1,E1,
     &	udf,vdf,pdf,cv_ksi,cv_eta, x, y, 
     &  yacobi, x_ksi, x_eta, y_ksi, y_eta )
	
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)  
	common /A1/U0,D,TAU 
        common /A6/uresid,vresid,presid
	common /A2/iXmax,iXsmax,iYmax,iYsmax
	common /timer/start, stoptime

        dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)
	dimension udf(iXm,iYm),vdf(iXm,iYm),pdf(iXm,iYm)
	dimension w(0:8),ex(0:8),ey(0:8)
	dimension Palph(-1:iXm+2,-1:iYm+2,0:8)
        dimension Ealph(0:iXm+1,0:iYm+1,0:8)

	dimension P1(-1:iXm+2,-1:iYm+2,0:8)
	dimension E1(0:iXm+1,0:iYm+1,0:8)

!     Grid data
        dimension x(iXm,iYm)
        dimension y(iXm,iYm)
	dimension yacobi(iXm,iYm)
        dimension x_ksi(iXm,iYm)
	dimension x_eta(iXm,iYm)
        dimension y_ksi(iXm,iYm)
	dimension y_eta(iXm,iYm)

!     Contravariant data array
	dimension cv_ksi(iXm,iYm,0:8)
	dimension cv_eta(iXm,iYm,0:8)

	open(1,file="fdlbm-2D-cy-work-data.csv")
	rewind(1)
        open(2,file="Cd-Cl-2RK2UP-steady.plt")
        rewind(2)
        write(2,*)"Variables=time,Cd,Cl"

!	open(9,file="Vort-2RK2UP-steady.plt")
!        rewind(9)
      write(*,*) 'Maximum timestep?'
      read(*,*) L

        call cpu_time(start)
      do 1000 num=1, L

	call ACP(num,p,u,v,ex,ey,w,Palph,Ealph,P1,E1,cv_ksi,cv_eta)	
	call BC(u,v,p,w,Palph,Ealph,E1)
	call MS(u,v,ex,ey,Palph)
		
      call SH(num,u,v,p,udf,vdf,pdf)
       if(mod(num,1000) .eq. 0) then
          call cpu_time(stoptime)
          write(*,*)"1000 step :",stoptime-start,"secs"
          start =stoptime
          endif

	if(mod(num,1000) .eq. 0) then
        call CdCl(num,u,v,p,x,y,yacobi,x_ksi,x_eta,y_ksi,y_eta)
	endif
        

       	if(mod(num,5000) .eq. 0) then
        call MSO(p,u,v,x,y)
	endif

        total_resid = 0.5*(uresid + vresid)
	if( total_resid  .lt.  10.0**-6.0) then
	write(*,*)"num:",num 
	go to 1100
	else
	endif

 1000 continue

 1100 continue

        call CdCl(num,u,v,p,x,y,yacobi,x_ksi,x_eta,y_ksi,y_eta)
        call MSO(p,u,v,x,y)

	close(1)
        close(2)
        close(9)

      return
	end	
********************************************************
*    ADVECTION AND COLLISION PROCESS                   *
********************************************************
      subroutine ACP(num,p,u,v,ex,ey,w,Palph,Ealph,P1,E1,cv_ksi,cv_eta)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)  
	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax 
	common /A3/TAUINV,ROUINV,ROUCs2,ROUCs2INV
	common /A4/ dt

	dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)
	dimension Palph(-1:iXm+2,-1:iYm+2,0:8)
	dimension Ealph(0:iXm+1,0:iYm+1,0:8)

	dimension P1(-1:iXm+2,-1:iYm+2,0:8)
	dimension E1(0:iXm+1,0:iYm+1,0:8)

	dimension w(0:8),ex(0:8),ey(0:8)

!     Contravariant data array
	dimension cv_ksi(iXm,iYm,0:8)
	dimension cv_eta(iXm,iYm,0:8)

        ad = 1.0
	do k=0,8
	do j=1, iYmax
	do i=2, iXmax-3
	isgksi =  sign(ad, cv_ksi(i,j,k) ) 
	isgeta =  sign(ad, cv_eta(i,j,k) )

	sgksi =  sign(ad, cv_ksi(i,j,k) ) 
	sgeta =  sign(ad, cv_eta(i,j,k) )

	iup = i-isgksi
	jup = j-isgeta
	iup2= i-2*isgksi
	jup2= j-2*isgeta
	iup3=i+isgksi
	jup3=j+isgeta
!
	Palph(i,0,k)=Palph(i,iYmax-1,k)
	Palph(i,-1,k)=Palph(i,iYmax-2,k)
	Palph(i,iYmax+1,k)=Palph(i,2,k)
	Palph(i,iYmax+2,k)=Palph(i,3,k)

!	Virtual points for i=1
	Palph(0,j,k)=2.0*Palph(1,j,k)-Palph(2,j,k)
	Palph(-1,j,k)=2.0*Palph(0,j,k)-Palph(1,j,k)
!
	Ealph(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

!	1./6. = 0.1666667
!
!3rd order upwind
	P1(i,j,k) = Palph(i,j,k)
     &-0.5*dt*TAUINV*(Palph(i,j,k)-Ealph(i,j,k))
     &-0.5*dt*cv_ksi(i,j,k)*sgksi*0.1666667*(3.0*Palph(i,j,k)
     &-6.0*Palph(iup,j,k)+Palph(iup2,j,k)+2.*Palph(iup3,j,k))      
     &-0.5*dt*cv_eta(i,j,k)*sgeta*0.1666667*(3.0*Palph(i,j,k)
     &-6.0*Palph(i,jup,k)+Palph(i,jup2,k)+2.*Palph(i,jup3,k))
!
!2rd order upwind
!	P1(i,j,k) = Palph(i,j,k)
!     &-0.5*dt*TAUINV*(Palph(i,j,k)-Ealph(i,j,k))
!     &-0.5*dt*cv_ksi(i,j,k)*sgksi*0.5*(3.0*Palph(i,j,k)
!     &-4.0*Palph(iup,j,k)+Palph(iup2,j,k))      
!     &-0.5*dt*cv_eta(i,j,k)*sgeta*0.5*(3.0*Palph(i,j,k)
!     &-4.0*Palph(i,jup,k)+Palph(i,jup2,k))

!	Cut-off line BC
	P1(i,1,k)=P1(i,iYmax,k)

	end do
	end do
	end do
!
	do k=0,8
	do j=1, iYmax
!	Equilibrium BC
	i=1
	u(i,j)=0.0
	v(i,j)=0.0
	p(i,j)=p(i+1,j) 
	P1(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

	P1(i,j,k)=P1(i,j,k)+(Palph(i+1,j,k)-Ealph(i+1,j,k)) !Based on T.S. Zhao

	end do
	end do

	do k=0,8
	do j=0.25*iYmax, 0.75*iYmax
	i=iXmax 
	u(i,j)=U0
	v(i,j)=0.0
	p(i,j)=2.0*p(i-1,j)-p(i-2,j)
	P1(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))
	end do
	end do

!
	do k=0, 8
	do j=1, 0.25*iYmax-1
	i=iXmax
	u(i,j)=u(i-1,j)
	v(i,j)=0.0
	p(i,j)=ROUCs2
	P1(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

	end do
	end do
!
	do k=0,8
	do j=0.75*iYmax+1, iYmax-1
	i=iXmax
	u(i,j)=u(i-1,j)
	v(i,j)=0.0
	p(i,j)=ROUCs2
	P1(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))
 
	end do
	end do

!     Second Step in the RK intergration
        ad = 1.0
	do k=0,8
	do j=1, iYmax
	do i=2, iXmax-3
	isgksi =  sign(ad , cv_ksi(i,j,k) ) 
	isgeta =  sign(ad , cv_eta(i,j,k) )

	sgksi =  sign(ad , cv_ksi(i,j,k) ) 
	sgeta =  sign(ad , cv_eta(i,j,k) )

	iup = i-isgksi
	jup = j-isgeta
	iup2= i-2*isgksi
	jup2= j-2*isgeta
	iup3=i+isgksi
	jup3=j+isgeta
!
	P1(i,0,k)=P1(i,iYmax-1,k)
	P1(i,-1,k)=P1(i,iYmax-2,k)
	P1(i,iYmax+1,k)=P1(i,2,k)
	P1(i,iYmax+2,k)=P1(i,3,k)

!	Virtual points for i=1
	P1(0,j,k)=2.0*P1(1,j,k)-P1(2,j,k)
	P1(-1,j,k)=2.0*P1(0,j,k)-P1(1,j,k)
!
!     Equilibrium function at mid-step; interpolated instead of updated
	E1(i,j,k) = 2.*Ealph(i,j,k) - E1(i,j,k) !odd? should be 2.0*E1(i,j,k)-Ealph(i,j,k)


!	1./6. = 0.1666667
!
!3rd Upwind
	Palph(i,j,k) = Palph(i,j,k)
     &-dt*TAUINV*(P1(i,j,k)-E1(i,j,k))
     &-dt*cv_ksi(i,j,k)*sgksi*0.1666667*(3.0*P1(i,j,k)
     &-6.0*P1(iup,j,k)+P1(iup2,j,k)+2.*P1(iup3,j,k))      
     &-dt*cv_eta(i,j,k)*sgeta*0.1666667*(3.0*P1(i,j,k)
     &-6.0*P1(i,jup,k)+P1(i,jup2,k)+2.*P1(i,jup3,k))

!2rd order upwind
!	Palph(i,j,k) = Palph(i,j,k)
!     &-dt*TAUINV*(P1(i,j,k)-E1(i,j,k))
!     &-dt*cv_ksi(i,j,k)*sgksi*0.5*(3.0*P1(i,j,k)
!     &-4.0*P1(iup,j,k)+P1(iup2,j,k))      
!     &-dt*cv_eta(i,j,k)*sgeta*0.5*(3.0*P1(i,j,k)
!     &-4.0*P1(i,jup,k)+P1(i,jup2,k))


!	Cut-off line BC
	Palph(i,1,k)=Palph(i,iYmax,k)

	end do
	end do
	end do

      return
	end

********************************************************
*     BOUNDARY CONDITION                               *
********************************************************
	subroutine BC(u,v,p,w,Palph,Ealph,E1)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)  

	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax
	common /A3/TAUINV,ROUINV,ROUCs2,ROUCs2INV

	dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)
	dimension Palph(-1:iXm+2,-1:iYm+2,0:8)
	dimension Ealph(0:iXm+1,0:iYm+1,0:8)
	dimension E1(0:iXm+1,0:iYm+1,0:8)
	dimension w(0:8),ex(0:8),ey(0:8)
        
        do k=0,8
        do j=1, iYmax
        do i=1,iXmax
        Ealph(i,j,k) = E1(i,j,k)   
        end do
        end do
        end do

	do k=0,8
	do j=1, iYmax
!	Equilibrium BC
	i=1
	u(i,j)=0.0
	v(i,j)=0.0
	p(i,j)= p(i+1,j) !Zero pressure gradient condition 
	Palph(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

	Palph(i,j,k)=Palph(i,j,k)+(Palph(i+1,j,k)-Ealph(i+1,j,k)) !Based on T.S. Zhao

	end do
	end do

	do k=0,8
	do j=0.25*iYmax, 0.75*iYmax
	i=iXmax 
	u(i,j)=U0
	v(i,j)=0.0
	p(i,j)=2.0*p(i-1,j)-p(i-2,j)
	Palph(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

	end do
	end do

	do k=0, 8
	do j=1, 0.25*iYmax-1
	i=iXmax
	u(i,j)=u(i-1,j)
	v(i,j)=0.0
	p(i,j)=ROUCs2
	Palph(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

	end do
	end do

	do k=0,8
	do j=0.75*iYmax+1, iYmax-1
	i=iXmax
	u(i,j)=u(i-1,j)
	v(i,j)=0.0
	p(i,j)=ROUCs2
	Palph(i,j,k)=w(k)*(p(i,j)+D*((ex(k)*u(i,j)+ey(k)*v(i,j))
     &              +1.5d0*(ex(k)*u(i,j)+ey(k)*v(i,j))**2.0d0
     &              -0.5d0*(u(i,j)**2.0d0+v(i,j)**2.0d0)))

	end do
	end do

!	Internal pressure
      do 320 j=1,iYmax
	  do 320 i=1,iXmax
	  p(i,j)=0.0d0
	    do 320 k=0,8
      p(i,j)=p(i,j)+Palph(i,j,k)
  320 continue

      return
	end

*******************************************
*     macroscopic status                  *
*******************************************
      subroutine MS(u,v,ex,ey,Palph)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)
	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax
	common /A3/TAUINV,ROUINV,ROUCs2,ROUCs2INV

	dimension Palph(-1:iXm+2,-1:iYm+2,0:8)
	dimension w(0:8),ex(0:8),ey(0:8)
	dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)

	do 410 j=1,iYmax
	  do 410 i=2,iXmax
	    u(i,j)=0.0d0
	    v(i,j)=0.0d0
	    do 400 k=0,8
		  u(i,j)=u(i,j) + ROUCs2INV*ex(k)*Palph(i,j,k)
		  v(i,j)=v(i,j) + ROUCs2INV*ey(k)*Palph(i,j,k)
  400 continue

  410 continue
      
      return
	end
****************************************
*     shuusoku hantei                  *
****************************************
      subroutine SH(num,u,v,p,udf,vdf,pdf)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300)
	common /A2/iXmax,iXsmax,iYmax,iYsmax
	common /A6/uresid,vresid,presid
        common /A7/uresid2,vresid2,presid2

	dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)
	dimension udf(iXm,iYm),vdf(iXm,iYm),pdf(iXm,iYm)

!     Local dimension
	dimension ua(iXmax,iYmax)
        dimension va(iXmax,iYmax)
        dimension pa(iXmax,iYmax)
!     
        usum2=0.0
	vsum2=0.0
	psum2=0.0
	iDomain = iXmax*0.25
        bunbo=dfloat(iDomain*iYmax)
!
	if(num.eq.1) then
	uresid=100.0d0
	vresid=100.0d0
	presid=100.0d0
	go to 600
	elseif(num.eq.2) then
      
	do 610 j=1,iYmax
	  do 610 i=1,iDomain
        udf(i,j)=(u(i,j)-ua(i,j))**2
	  vdf(i,j)=(v(i,j)-va(i,j))**2
	  pdf(i,j)=(p(i,j)-pa(i,j))**2
        usum2=usum2+udf(i,j)
	  vsum2=vsum2+vdf(i,j)
	  psum2=psum2+pdf(i,j)

  610 continue
      uresid2=dsqrt(usum2/bunbo)
	vresid2=dsqrt(vsum2/bunbo)
	presid2=dsqrt(psum2/bunbo)
	uresid=uresid2/uresid2
	vresid=vresid2/vresid2
	presid=presid2/presid2
	write(*,*)num,uresid,vresid,presid
	else
      do 620 j=1,iYmax
	  do 620 i=1,iDomain
        udf(i,j)=(u(i,j)-ua(i,j))**2
	  vdf(i,j)=(v(i,j)-va(i,j))**2
	  pdf(i,j)=(p(i,j)-pa(i,j))**2
        usum2=usum2+udf(i,j)
	  vsum2=vsum2+vdf(i,j)
	  psum2=psum2+pdf(i,j)
  620 continue
        uresid=sqrt(usum2/bunbo)/uresid2
	vresid=sqrt(vsum2/bunbo)/vresid2
	presid=sqrt(psum2/bunbo)/presid2
	endif

	do 630 num1=1,10000
	if(num.eq.num1*1000) then

	write(*,*) num, uresid, vresid, presid
!	write(100,*) num, uresid, vresid, presid

	endif
  630 continue

!	if(mod(num,1000) .eq. 0) then
!	write(*,*) num, uresid, vresid, presid
!	write(100,*) num, uresid, vresid, presid
!	endif
 
  600 continue
      do 640 j=1,iYmax
	  do 640 i=1,iDomain
      ua(i,j)=u(i,j)
	va(i,j)=v(i,j)
	pa(i,j)=p(i,j)
  640 continue
      return
	end

*******************************************
*     macro status output                 *
*******************************************
      subroutine MSO(p,u,v,x,y)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300) 
	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax 
	dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)

!Grid data array     
	dimension x(iXm,iYm)
        dimension y(iXm,iYm)

!     Local dimension
        dimension pstd(iXm,iYm) 


!
        pmean = 0.0
        do j=1, iYmax
        do i=1, iXmax
        pmean = pmean + p(i,j)   
        end do
        end do   
        pmean = pmean/(iXmax*iYmax)

        do j=1, iYmax
        do i=1, iXmax
        pstd(i,j) = (p(i,j) - pmean)/pmean
        pstd(i,j) = pstd(i,j)*1000.0
        end do
        end do
!
!
!	write(1,*) "Variables = x, y, p, u, v,pstd"
!	write(1,*) "Zone i=",iXmax-3, "j=",iYmax
!	do 500 j=1,iYmax
!	  do 500 i=1,iXmax-3
!        write(1,*)  x(i,j), y(i,j), p(i,j),u(i,j),v(i,j), pstd(i,j)
! 500  continue

!     writing output for paraview in csv format
 900    format (F12.4, A, F12.4, A, F12.4, A, F12.4, A, F12.4, A, F12.4)
 902    format (F12.4, A, F12.4, A, F12.4)
 901    format (6E12.4)
        bb = 0.0
        write(1,*) "x,y,z,p,u,v"
!        do j = 1, iYmax
         do i = 1, iXmax-3
         do j = 1, iYmax
            write(1, 900)  x(i,j),",",y(i,j),",",bb,",",u(i,j),",",
     & v(i,j),",",p(i,j)
        end do
        end do
 

	return
        end

!***************************************************************
       subroutine CdCl(num,u,v,p,x,y,yacobi,x_ksi,x_eta,y_ksi,y_eta)
	implicit real*8(a-h,o-z)
	parameter(iXm=300,iYm=300) 

	common /A1/U0,D,TAU
	common /A2/iXmax,iXsmax,iYmax,iYsmax 
	common /A4/ dt

!Local common variables
!        common /angle/ ds
!        common /qty/ Cd, Cl, Fx, Fy

	dimension u(iXm,iYm),v(iXm,iYm),p(iXm,iYm)

!     Grid data array
	dimension x(iXm,iYm)
        dimension y(iXm,iYm)
	dimension yacobi(iXm,iYm)
        dimension x_ksi(iXm,iYm) 
	dimension x_eta(iXm,iYm) 
        dimension y_ksi(iXm,iYm)
	dimension y_eta(iXm,iYm)

!     Local dimension
        dimension uzu(iXmax,iYmax)
        dimension us(0:iXmax+1,0:iYmax+1)
        dimension vs(0:iXmax+1,0:iYmax+1)

        dimension shear(1,iYmax) !For seperation angle

        Fx=0.0
        Fy=0.0
        Cd=0.0
        Cl=0.0

        do j=1,iYmax-1
        i=1   
        Fx = Fx - 0.5*(p(i,j)+p(i,j+1))*x(i,j) 
     & *sqrt( (( x(i,j+1)-x(i,j) ))**2 + (( y(i,j+1)-y(i,j) ))**2)   
!Inserting viscosity effect
        Fx = Fx + 2.0*D*(TAU/3.0)*y(i,j) 
     & *sqrt( (( x(i,j+1)-x(i,j) ))**2 + (( y(i,j+1)-y(i,j) ))**2) 
     &*0.5*((-x_eta(i,j)/yacobi(i,j))*(4.0*u(i+1,j)-3.0*u(i,j)-u(i+2,j))
     & +(x_ksi(i,j)/yacobi(i,j))*(u(i,j+1)-u(i,j-1)) 
     &+ (y_eta(i,j)/yacobi(i,j))*(4.0*v(i+1,j)-3.0*v(i,j)-v(i+2,j)) 
     &+(-y_ksi(i,j)/yacobi(i,j))*(v(i,j+1)-v(i,j-1)) )
        end do

       do j=1,iYmax-1
        i=1   
        Fy = Fy - 0.5*(p(i,j)+p(i,j+1))*y(i,j) 
     & *sqrt( (( x(i,j+1)-x(i,j) ))**2 + (( y(i,j+1)-y(i,j) ))**2)   

        Fy = Fy + 2.0*D*(TAU/3.0)*x(i,j) 
     & *sqrt( (( x(i,j+1)-x(i,j) ))**2 + (( y(i,j+1)-y(i,j) ))**2) 
     &*0.5*((-x_eta(i,j)/yacobi(i,j))*(4.0*u(i+1,j)-3.0*u(i,j)-u(i+2,j))
     & +(x_ksi(i,j)/yacobi(i,j))*(u(i,j+1)-u(i,j-1)) 
     &+ (y_eta(i,j)/yacobi(i,j))*(4.0*v(i+1,j)-3.0*v(i,j)-v(i+2,j)) 
     &+(-y_ksi(i,j)/yacobi(i,j))*(v(i,j+1)-v(i,j-1)) )
        end do



      Cd=abs(Fx)/(D*U0**2*1.0) !radius=1.0  
      Cl=Fy/(D*U0**2*1.0) !radius=1.0  

      ad = 1.0
      
      do j=2, iYmax/4-1
      i=2
      if( sign(ad, v(i,j) )  >= 0.0 ) then
      shear(1,j) = y(i,j)   
      endif
      end do

      ds = maxval(shear, mask = shear <=1.0)
      ds = (180*asin(ds))/3.14159

!      write(*,*)"Fx  : ", Fx
!      write(*,*)"Fy  : ", Fy
      write(*,*) "Cd : ", Cd
      write(*,*) "Cl : ", Cl
      write(*,*) "num : ", num

      write(2,*) num*dt, Cd, Cl

!      if(mod(num,5000) .eq. 0)then
!     Calculate vorticity field
!     Initiating array value
      do j=1,iYmax
      do i=1,iXmax
         us(i,j) = u(i,j)
         vs(i,j) = v(i,j)
      end do
      end do

!	Vorticity in general grid calculation
!	do j=1, iYmax-1
!	do i=2, iXmax-3
!	us(i,0) = u(i,iYmax)
!	vs(i,0) = v(i,iYmax)
!	uzu(i,j) = ( ((y_eta(i,j)/yacobi(i,j))*0.5*(vs(i,j)
!     &-vs(i-1,j))))
!     &+(-(y_ksi(i,j)/yacobi(i,j))*0.5*(vs(i,j+1)-vs(i,j-1)))
!     &-((-(x_eta(i,j)/yacobi(i,j))*0.5*(us(i+1,j)-us(i-1,j)))
!     &+((x_ksi(i,j)/yacobi(i,j))*0.5*(us(i,j+1)-us(i,j-1))) )

!        uzu(i,iYmax)=uzu(i,1)

!	end do
!	end do


!	do j=1,iYmax-1
!	i=1
!	dvdksi=	-0.5*(3.0d0*vs(i,j)-4.0d0*vs(i+1,j)+vs(i+2,j))
!	dvdeta=0.5*(vs(i,j+1)-vs(i,j-1))

!	dudksi= -0.5*(3.0d0*us(i,j)-4.0d0*us(i+1,j)+us(i+2,j))
!	dudeta =0.5*(us(i,j+1)-us(i,j-1))
!	uzu(i,j)= ( y_eta(i,j)/yacobi(i,j) )*dvdksi
!     &	+( -y_ksi(i,j)/yacobi(i,j) )*dvdeta
!     &-( (-x_eta(i,j)/yacobi(i,j))*dudksi
!     & +(x_ksi(i,j)/yacobi(i,j))*dudeta ) 

!        uzu(i,iYmax)=uzu(i,1)

!	end do
        

!	write(9,*) "Variables = x, y,u,v,vort"
!	write(9,*) "Zone i=", iXmax-3, "j=", iYmax
!	do 10 j=1,iYmax
!	do 10 i=1,iXmax-3 
!	write(9,*) x(i,j), y(i,j),u(i,j),v(i,j), uzu(i,j)
!   10 continue

!      endif

        
        return
        end
