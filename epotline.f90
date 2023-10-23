program epotline
  implicit none !Auto to kanoume panta oste na kanoume declaration of variables oste na min orizontai mesa sto programma aftereta (implicit kanones tis fortran diaforoi me tis python)
  !--------------------------------------------------------------------------------------------
  !Declaration of variables:
  real, allocatable     :: X(:),Y(:),Q(:)    !Dimiourgoume ta array mas, opou gia poio sosto memory allocation tha vroume tin diastasi tous pio kato.
  integer               :: N                 !O arithmos ton fortion mesa sto provlima
  real                  :: x0,y0             !Simeio ypolgismou P(x,y)
  integer               :: ndraw,i,j         !The ndraw is for the number of the lines to draw and i,j is for do loops.
  real                  :: theta             !H gonia ypologismou tis kyklikis epifaneias gyro apo to fortio
  real,parameter        :: PI=3.141592654    !Exoume san parameter (static variable) tin stathera pi gia tous trigonometrikous ypologismous
  real                  :: L                 !Exoume mia metavliti typou real pou mas orizei to tetragono sxediasmou ton isodynamikon epifaneion
  real                  :: rmax,rmin         !Gia na kseroume kathe fora poso konta kai makria eimaste sto pio kontino kai makrino fortio antistoixa
  !--------------------------------------------------------------------------------------------
  !User Interface:
  print *, '#Enter the number of charges:'
  read  *, N
  print *, '#N= ',N
  !Array Allocation:
  ALLOCATE(X(N))
  ALLOCATE(Y(N))
  ALLOCATE(Q(N))
  !----------------------------------------
  do i=1,N,1
     print *, '#Charge No: ',i
     print *, '#Give the position (x,y) and the charge for this situation: '
     read  *, X(i),Y(i),Q(i)
     print *, '#(x,y)= ',X(i),Y(i), '#Q= ',Q(i)
  end do
  !----------------------------------------
  !Calculations:
  !We draw the epotlines in the square -L<=x<=+L and -L<=y<=+L
  ndraw=4 !We give an initial value
  L=1.0   !Also an initial value
  do i=-ndraw,ndraw,1
     do j=-ndraw,ndraw,1
        x0=i*(L/ndraw)
        y0=j*(L/ndraw)
        call mdist(x0,y0,X,Y,N,rmin,rmax)
        !Apofeygoume na pame konta se kapoio fortio
        if (rmin >= L/(10*ndraw)) then
           call potencial(x0,y0,X,Y,Q,N)
        end if
        !---------------------------------------------------------------------
     end do
     !-------------------------------------------------------------------------
  end do
  !----------------------------------------------------------------------------
end program epotline
!===============================================================================
subroutine mdist(x0,y0,X,Y,N,rmin,rmax)
  implicit none
  !-----------------------------------------------------------------------------
  !Delcaration of variables:
  integer           :: N
  real,dimension(N) :: X,Y
  real              :: x0,y0
  real              :: rmin,rmax
  integer           :: i
  real              :: r  !Dummy index for every calculation of distance
  !-----------------------------------------------------------------------------
  !Calculations:
  rmin=HUGE(rmin)
  rmax=-HUGE(rmax) !The HUGE(A) function gives the largest number for the type of variables A
  do i=1,N,1
     r=sqrt((x0-X(i))**2 + (y0-Y(i))**2)
     if (r>rmax) then
        rmax=r
     end if
     !-------------------------------------------------------------------------
     if (r<rmin) then
        rmin=r
     end if
     !--------------------------------------------------------------------------
  end do
  !-----------------------------------------------------------------------------
end subroutine mdist
!================================================================================
subroutine efield(x0,y0,X,Y,Q,N,Ex,Ey)
  implicit none
  !------------------------------------------------------------------------------
  !Declaration of variables:
  integer              :: N
  real,dimension(N)    :: X,Y,Q
  real                 :: x0,y0
  real                 :: Ex,Ey
  integer              :: i
  real                 :: r3,xi,yi !r3 is the distance ((x0-x(i))^2+(y0-Y(i))^2)^(-3/2) and xi and yi are the inside difference between points
  !-----------------------------------------------------------------------------
  !Calculations:
  Ex=0.0
  Ey=0.0  !We give the initial value to the electric field before calculations
  do i=1,N,1
     xi=x0-X(i)
     yi=y0-Y(i)
     r3=(xi*xi+yi*yi)**(-1.5) !The multiplication is a better operation from division in fortran
     Ex=Ex+Q(i)*xi*r3
     Ey=Ey+Q(i)*yi*r3
  end do
  !------------------------------------------------------------------------------
end subroutine efield
!================================================================================
subroutine potencial(xin,yin,X,Y,Q,N)
  implicit none
  !------------------------------------------------------------------------------
  !Delcaration of variables:
  integer           :: N                 !The number of charges in our distirbution
  real,dimension(N) :: X,Y,Q             !Our arrays for our problem.Also we need to declare first the N integer to have the size of our arrays.
  real              :: xin,yin           !Our initial values that show the point P(x,y) we calculate the dynamic line
  real              :: x0,y0             !Local variables that expresses all the points of our dynamic lines
  real,parameter    :: step=0.01         !It's the dl step every time (given from the teacher a good static value for it)
  real,parameter    :: max_distance=20.0 !The maximum distance away from a charge
  integer           :: i,direction       !The i variable is for do loops and direction is an integer with values +1 or -1 that indicates our direction on the dynamic line
  real              :: rmin,rmax         !The minimum and maximum distance from our point in dynamic line to the charges
  real              :: r,dx,dy,dl        !The distance r from the starting point xin,yin every time we move in the lattice, dx and dy the steps and dl is known
  real              :: Ex,Ey,E           !The electric field values
  real              :: check1,check2     !It holds the previous values of dx and dy
  integer           :: count             !It counts the change of direction when we draw a single dynamic line
  !------------------------------------------------------------------------------
  !Calculations:
  !Initial Values:
  dl=step              !We have that dl takes the initial value of step because epotiles are closed curves
  x0=xin
  y0=yin               !Initial values given from the user
  dx=0.0
  dy=0.0               !Initial values for dx and dy, before calculations
  r=step               !We give an initial value in order to start a loop
  check1=dx
  check2=dy
  count=0
  do while (r>=(0.9*dl) .and. r<=max_distance)
     print *, x0,y0                   !We print the initial value of our current place
     call efield(x0+0.5*dx,y0+0.5*dy,X,Y,Q,N,Ex,Ey) !We calculate Ex and Ey for our current point x0,y0.Akomi vazoume to x0+0.5*dx y0+0.5*dy oste na meiosoume to sfalma tis diakritopoiisis kai na exoume tin mesi syneisfora apo to ilektriko pedio.
     E=sqrt(Ex*Ex+Ey*Ey)
     !--------------------------------------------------------------------
     !Check when we have E=0 or very close to 0
     if (E<1.0e-10) then
        exit !Stop the current do loop and moves to the next do loop of direction
     end if
     !---------------------------------------------------------------------
     dx=dl*(Ey/E) !If E=0 we have a big problem.But it solved with the line of code above.
     dy=-dl*(Ex/E) !The calculations of our step
     !We check if we had a change of direction during our draw of the single dynamic line
     if (check1*dx<0.0 .or. check2*dy<0.0) then
        count=count+1
     end if
     !---------------------------------------------------------------------
     x0=x0+dx
     y0=y0+dy
     !We calculate now the distance of x0,y0 from the starting point xin,yin
     r=sqrt((x0-xin)**2 + (y0-yin)**2)
     check1=dx
     check2=dy
     !We solve the problem of putting equal charges to a polygone
     !If the nymber of changes of direction is more than 2*N then we fell into an infinite dynamic well (apeirovatho pigadi dynamikou)
     if (count>N**4) then
        exit
     end if
     !----------------------------------------------------------------------
  end do
  !------------------------------------------------------------------------------
end subroutine potencial
!================================================================================
