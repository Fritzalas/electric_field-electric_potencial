program elines
  implicit none !Auto to kanoume panta oste na kanoume declaration of variables oste na min orizontai mesa sto programma aftereta (implicit kanones tis fortran diaforoi me tis python)
  !--------------------------------------------------------------------------------------------
  !Declaration of variables:
  real, allocatable     :: X(:),Y(:),Q(:)    !Dimiourgoume ta array mas, opou gia poio sosto memory allocation tha vroume tin diastasi tous pio kato.
  integer               :: N                 !O arithmos ton fortion mesa sto provlima
  real                  :: x0,y0             !Simeio ypolgismou P(x,y)
  integer               :: ndraw,i,j         !The ndraw is for the number of the lines to draw and i,j is for do loops.
  real                  :: theta             !H gonia ypologismou tis kyklikis epifaneias gyro apo to fortio
  real,parameter        :: PI=3.141592654    !Exoume san parameter (static variable) tin stathera pi gia tous trigonometrikous ypologismous
  integer               :: find_nd           !Kanoume declare tin function mas
  real                  :: k                 !Gia tin xrisi tou sto function
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
  do i=1,N,1
     k=Q(i)
     ndraw=find_nd(k,Q,N)
     print *, '#The current ndraw is equal to:',ndraw
     do j=1,2*ndraw,1
        theta=(PI/ndraw)*j
        x0=X(i)+0.1*cos(theta)
        y0=Y(i)+0.1*sin(theta)
        call eline(x0,y0,X,Y,Q,N)
     end do
  end do
  !-----------------------------------------
end program elines
!Vazoume ena sxolio oste o metaglotistis na diavasei pio kato apo to end program pou kapoies fores to paraleipei kai vgazei error
!=======================================================================================
function find_nd(k,Q,N)
  implicit none
  !--------------------------------------------------
  !Declaration of variables:
  integer            :: N
  real,dimension(N)  :: Q
  real               :: k          !The local charge we visit every time
  integer            :: find_nd    !The result of our function gives the numbers of dynamic lines drawn for each charge
  real               :: min_charge !The minimum charge for our given array
  integer            :: i
  real               :: current_absolute_value
  !------------------------------------------------------------------------------
  !Calculations:
  !Initialize min_absolute_value with the absolute value of the first element of the array
  min_charge=abs(Q(1))
  ! Iterate through the array to find the minimum absolute value
  do i =2,N,1
    current_absolute_value = abs(Q(i))
    if (current_absolute_value < min_charge) then
      min_charge = current_absolute_value
    end if
  end do
  find_nd=NINT(abs(k)/min_charge)*6
end function find_nd
!================================================================================
subroutine eline(xin,yin,X,Y,Q,N)
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
  real              :: r,dx,dy,dl        !The distance r from every charge, dx and dy the steps and dl is known
  real              :: Ex,Ey,E           !The electric field values
  real              :: check1,check2     !It holds the previous values of dx and dy
  integer           :: count             !It counts the change of direction when we draw a single dynamic line
  !------------------------------------------------------------------------------
  !Calculations:
  !Both Directions:
  do direction=-1,1,2     !Step here is 2 so we get -1 and +1 for both directions
     dl=step*direction    !We have that dl goes to both directions. When dl is positive we move with the vectors of electric field and when it's negative to the opposite direction.
     x0=xin
     y0=yin               !Initial values given from the user
     dx=0.0
     dy=0.0               !Initial values for dx and dy, before calculations
     check1=dx
     check2=dy
     count=0              !Initial values for our new variables
     !Now we want the maximum and minimum distance from every charge
     call mdist(x0,y0,X,Y,N,rmin,rmax) !mdist subroutine gives us the rmin and rmax every time for the given x0,y0
     !Check about space limits for calculations:
     do while (rmin>2.0*step .and. rmax<max_distance)
        print *, x0,y0                   !We print the initial value of our current place
        call efield(x0+0.5*dx,y0+0.5*dy,X,Y,Q,N,Ex,Ey) !We calculate Ex and Ey for our current point x0,y0.Akomi vazoume to x0+0.5*dx y0+0.5*dy oste na meiosoume to sfalma tis diakritopoiisis kai na exoume tin mesi syneisfora apo to ilektriko pedio.
        E=sqrt(Ex*Ex+Ey*Ey)
        !--------------------------------------------------------------------
        !Check when we have E=0 or very close to 0
        if (E<1.0e-10) then
           exit !Stop the current do loop and moves to the next do loop of direction
        end if
        !---------------------------------------------------------------------
        dx=dl*(Ex/E) !If E=0 we have a big problem.But it solved with the line of code above.
        dy=dl*(Ey/E) !The calculations of our step
        !We check if we had a change of direction during our draw of the single dynamic line
        if (check1*dx<0.0 .or. check2*dy<0.0) then
           count=count+1
        end if
        !---------------------------------------------------------------------
        x0=x0+dx 
        y0=y0+dy     !We move dx and dy as we saw in our theory
        check1=dx
        check2=dy
        !We solve the problem of putting equal charges to a polygone
        !If the nymber of changes of direction is more than 2*N then we fell into an infinite dynamic well (apeirovatho pigadi dynamikou)
        if (count>2*N) then
           exit
        end if
        !----------------------------------------------------------------------
        call mdist(x0,y0,X,Y,N,rmin,rmax) !We call mdist to make again for the next loop the neccessary check
     end do
     !-------------------------------------------------------------------------
  end do
  !-----------------------------------------------------------------------------
end subroutine eline
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
