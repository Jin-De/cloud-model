PROGRAM dry
IMPLICIT NONE

INTEGER, PARAMETER :: nx=83, nz=42
REAL, PARAMETER :: g=9.8, Cp=1004, Cv=717, Rd=287.047 ,p0=100000, ps=96500,  &
dx=400, dz=400, dt=2, dtx=2*dt/dx, dtz=2*dt/dz, delta=3, zcnt=3000, imid=43, &
radx=4000, radz=4000, tripi=4.*atan(1.0), Cs=50, kx=50, kz=50
INTEGER :: i, k, j
REAL, DIMENSION(nz) :: tb, qb, pb, pib, rhou, rhow
REAL, DIMENSION(nx,nz) :: thp, th, thm, dth,  &
                          up, u, um, du,      &
                          wp, w, wm, dw,      &
                          pip, pi, pim, dpi, df   
REAL :: rad, tave, z_scalar, tup, tdn

! mean state initialize
tb(2)=300.
qb(2)=0.
pb(2)=ps
pib(2)=(ps/p0)**(Rd/Cp)
rhou(2)=(p0*(pib(2))**(Cv/Rd))/(Rd*tb(2))
rhow(2)=rhou(2)
DO k=3,nz-1
  tb(k)=300.
  qb(k)=0.
  tave=0.5*(tb(k)+tb(k-1))
  pib(k)=pib(k-1)-g*dz/(Cp*tave) 
  pb(k)=p0*(pib(k))**(Cp/Rd)
  rhou(k)=(p0*(pib(k))**(Cv/Rd))/(Rd*tave)
  rhow(k)=0.5*(rhou(k)+rhou(k-1))
ENDDO
tb(1)=tb(2)
tb(nz)=tb(nz-1)
qb(1)=qb(2)
qb(nz)=qb(nz-1)
pib(1)=pib(2)
pib(nz)=pib(nz-1)
pb(1)=pb(2)
pb(nz)=pb(nz-1)
rhou(1)=rhou(2)
rhou(nz)=rhou(nz-1)
rhow(1)=rhow(2)
rhow(nz)=rhow(nz-1)

th=0
u=0
um=0
w=0
wm=0
pi=0
! perturbation initialization
DO i=2,nx-1
  DO k=2,nz-1
    z_scalar=dz*(real(k)-0.5)
    rad=sqrt(((z_scalar-zcnt)/radz)**2+(dx*(real(i)-imid)/radx)**2)
    IF (rad>=1)THEN
      th(i,k)=0
    ELSE
      th(i,k)=0.5*delta*(cos(tripi*rad)+1)
    ENDIF
    thm(i,k)=th(i,k)
  ENDDO
ENDDO

DO i=2,nx-1
  pi(i,nz)=0
  DO k=nz-1,2,-1
    tup=th(i,k+1)/(tb(k+1)**2)
    tdn=th(i,k)/(tb(k)**2)
    pi(i,k)=pi(i,k+1)-0.5*(g/Cp)*(tup+tdn)*dz
    pim(i,k)=pi(i,k)
  ENDDO
  pi(i,1)=pi(i,2)
  pim(i,1)=pim(i,2)
ENDDO

OPEN(10,FILE='initial_condition.dat',access='direct',recl=nz*4*4)
WRITE(10,rec=1) tb, pib, pb, rhou
CLOSE(10)

OPEN(10,FILE='result.dat',access='direct',recl=(nx-2)*(nz-2)*4*4)

! integration
DO j=1,2000

  call u_equation(du,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)
  call diffusion(df,u,kx,kz,dx,dz,dt,nx,nz)
  up=um+du+df
  call w_equation(dw,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)
  call diffusion(df,w,kx,kz,dx,dz,dt,nx,nz)
  wp=wm+dw+df
  call th_equation(dth,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)
  call diffusion(df,th,kx,kz,dx,dz,dt,nx,nz)
  thp=thm+dth+df
  call pi_equation(dpi,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)  
  call diffusion(df,pi,kx,kz,dx,dz,dt,nx,nz)
  pip=pim+dpi+df

  ! boundary
  wp(:,nz-1)=0
  wp(:,2)=0
  call boundary(up,nx,nz)
  call boundary(wp,nx,nz)
  call boundary(thp,nx,nz)
  call boundary(pip,nx,nz)
  
   um =   u
    u =  up
   wm =   w
    w =  wp
  thm =  th
   th = thp
  pim =  pi
   pi = pip
  
  ! output
  WRITE(10,rec=j) u(2:nx-1,2:nz-1),w(2:nx-1,2:nz-1), &
                  th(2:nx-1,2:nz-1),pi(2:nx-1,2:nz-1)
ENDDO

CLOSE(10)

END PROGRAM dry

SUBROUTINE u_equation(du,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)

INTEGER :: i, k
INTEGER, INTENT(in) :: nx, nz
REAL, PARAMETER :: Cp=1004
REAL, INTENT(in) :: dtx, dtz, dt
REAL, INTENT(in), DIMENSION(nz) :: tb, rhou, rhow
REAL, INTENT(in), DIMENSION(nx,nz) :: u, w, th, pi
REAL, INTENT(out), DIMENSION(nx,nz) :: du

du=0
DO i=2,nx-1
  DO k=2,nz-1
    du(i,k)= -0.25 * dtx * ( (u(i+1,k) + u(i,k) )**2 - ( u(i,k) + u(i-1,k) )**2) &
             -0.25*dtz*( &
                rhow(k+1) * ( w(i,k+1) + w(i-1,k+1) )*( u(i,k+1) + u(i,k  ) ) &
               -rhow(k  ) * ( w(i,k  ) + w(i-1,k  ) )*( u(i,k  ) + u(i,k-1) ) &
                )/rhou(k) &
             -dtx * Cp * tb(k) * ( pi(i,k) - pi(i-1,k) )
  ENDDO
ENDDO

RETURN
END SUBROUTINE u_equation

SUBROUTINE w_equation(dw,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)

INTEGER :: i, k
INTEGER, INTENT(in) :: nx, nz
REAL, PARAMETER :: Cp=1004, g=9.8
REAL, INTENT(in) :: dtx, dtz, dt
REAL, INTENT(in), DIMENSION(nz) :: tb, rhou, rhow
REAL, INTENT(in), DIMENSION(nx,nz) :: u, w, th, pi
REAL, INTENT(out), DIMENSION(nx,nz) :: dw

dw=0
DO i=2,nx-1
  DO k=2,nz-1
    dw(i,k)= -0.25*dtx*( &
                   ( u(i+1,k) + u(i+1,k-1) ) * ( w(i+1,k) + w(i  ,k) ) &
                  -( u(i  ,k) + u(i  ,k-1) ) * ( w(i  ,k) + w(i-1,k) ) ) &
             -0.25*dtz*( &
                   rhou(k  ) * ( w(i,k+1) + w(i,k  ) )**2 &
                  -rhou(k-1) * ( w(i,k  ) + w(i,k-1) )**2 ) &
                  /rhow(k) &
             -0.5*dtz*Cp* ( tb(k) + tb(k-1) ) * ( pi(i,k) - pi(i,k-1) ) &
             +dt * g * ( th(i,k)/tb(k) + th(i,k-1)/tb(k-1) )
  ENDDO
ENDDO

RETURN
END SUBROUTINE w_equation

SUBROUTINE th_equation(dth,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)

INTEGER :: i, k
INTEGER, INTENT(in) :: nx, nz
REAL, PARAMETER :: Cp=1004
REAL, INTENT(in) :: dtx, dtz, dt
REAL, INTENT(in), DIMENSION(nz) :: tb, rhou, rhow
REAL, INTENT(in), DIMENSION(nx,nz) :: u, w, th, pi
REAL, INTENT(out), DIMENSION(nx,nz) :: dth

dth=0
DO i=2,nx-1
  DO k=2,nz-1
    dth(i,k)= -0.5*dtx*( &
                   u(i+1,k) * ( th(i+1,k) + th(i  ,k) ) &
                  -u(i  ,k) * ( th(i  ,k) + th(i-1,k) ) ) &
              -0.5*dtz*( &
                   rhow(k+1) * w(i,k+1) * ( th(i,k+1) + th(i,k  ) ) &
                  -rhow(k  ) * w(i,k  ) * ( th(i,k  ) + th(i,k-1) )) &
                  /rhou(k) &
              -0.5*dtz*( &
                   rhow(k+1) * w(i,k+1) * ( tb(k+1) - tb(k  ) )  &
                  +rhow(k  ) * w(i,k  ) * ( tb(k  ) - tb(k-1) ))
  ENDDO
ENDDO

RETURN
END SUBROUTINE th_equation

SUBROUTINE pi_equation(dpi,u,w,th,pi,tb,rhou,rhow,dtx,dtz,dt,nx,nz)

INTEGER :: i, k
INTEGER, INTENT(in) :: nx, nz
REAL, PARAMETER :: Cp=1004, Cs=50
REAL, INTENT(in) :: dtx, dtz, dt
REAL, INTENT(in), DIMENSION(nz) :: tb, rhou, rhow
REAL, INTENT(in), DIMENSION(nx,nz) :: u, w, th, pi
REAL, INTENT(out), DIMENSION(nx,nz) :: dpi

dpi=0
DO i=2,nx-1
  DO k=2,nz-1
    dpi(i,k)= -Cs**2/(rhou(k)*Cp*tb(k)**2)*( &
                  dtx*rhou(k)*tb(k)* ( u(i+1,k) - u(i,k) ) &
                +0.5*dtz*( &
                   rhow(k+1) * w(i,k+1) * ( tb(k+1) + tb(k  ) ) &
                  -rhow(k  ) * w(i,k  ) * ( tb(k  ) + tb(k-1) ))) 
  ENDDO
ENDDO

RETURN
END SUBROUTINE pi_equation

SUBROUTINE boundary(var,nx,nz)

INTEGER, INTENT(in) :: nx, nz
REAL, INTENT(inout), DIMENSION(nx,nz) :: var

var(:,1)=var(:,2)
var(:,nz)=var(:,nz-1)
var(nx,:)=var(2,:)
var(1,:)=var(nx-1,:)

RETURN
END SUBROUTINE boundary

SUBROUTINE diffusion(dvar,var,kx,kz,dx,dz,dt,nx,nz)

INTEGER, INTENT(in) :: nx, nz
REAL, INTENT(in) :: dx, dz, dt, kx, kz
REAL, INTENT(in), DIMENSION(nx,nz) :: var
REAL, INTENT(out), DIMENSION(nx,nz) :: dvar

dvar=0
DO i=2,nx-1
  DO k=2,nz-1
    dvar(i,k)= 2*dt*( &
                 kx*(var(i+1,k  )-2*var(i,k)+var(i-1,k  ))/(dx**2) &
                +kz*(var(i  ,k+1)-2*var(i,k)+var(i  ,k-1))/(dz**2))
  ENDDO
ENDDO

RETURN
END SUBROUTINE diffusion

