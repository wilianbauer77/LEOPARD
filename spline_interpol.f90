!> Uses natural cubic splines to interpolate the array y over the grid x (routine was provided by Alex Godunov, Old Dominion University, Norfolk, VA, USA)
!! \param x array of data abscissas (in strictly increasing order)
!! \param y array of data ordinates
!! \param b array of spline coefficients
!! \param c array of spline coefficients
!! \param d array of spline coefficients
!! \param n size of the arrays x and y (n>=2)
subroutine spline_interpol (x, y, b, c, d, n)
  !http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
  !======================================================================
  !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
  !  for cubic spline interpolation
  !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !  for  x(i) <= x <= x(i+1)
  !  Alex G: January 2010
  !======================================================================
  implicit none
  integer :: n
  real, dimension(n) ::  x
  real, dimension(n) :: y, b, c, d
  integer :: i, j, gap
  real :: h

  gap = n-1

  ! check input
  if ( n < 2 ) return
  if ( n < 3 ) then
     b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
     c(1) = 0.0
     d(1) = 0.0
     b(2) = b(1)
     c(2) = 0.0
     d(2) = 0.0
     return
  end if

  !
  ! step 1: preparation
  !

  d(1) = x(2) - x(1)
  c(2) = (y(2) - y(1))/d(1)
  do i = 2, gap
     d(i) = x(i+1) - x(i)
     b(i) = 2 *(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
  end do

  !
  ! step 2: end conditions 
  !

  b(1) = -d(1)
  b(n) = -d(n-1)
  c(1) = 0.0
  c(n) = 0.0
  if(n /= 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n) = -c(n)*d(n-1)**2 /(x(n)-x(n-3))
  end if

  !
  ! step 3: forward elimination 
  !

  do i = 2, n
     h = d(i-1)/b(i-1)
     b(i) = b(i) - h*d(i-1)
     c(i) = c(i) - h*c(i-1)
  end do

  !
  ! step 4: back substitution
  !

  c(n) = c(n)/b(n)
  do j = 1, gap
     i = n-j
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
  end do

  !
  ! step 5: compute spline coefficients
  !

  b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2*c(n))
  do i = 1, gap
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3*c(i)
  end do

  c(n) = 3 *c(n)
  d(n) = d(n-1)

  do i=1,n
     if(y(i)==0.0) then
        b(i)=0.0
        c(i)=0.0
        d(i)=0.0
     endif

  enddo

end subroutine spline_interpol
