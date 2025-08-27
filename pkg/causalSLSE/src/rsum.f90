subroutine exprsum(x, n, value)
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)    :: n
  real(c_double), intent(in)    :: x(n)
  real(c_double), intent(out)   :: value

  integer                      :: i
  integer, parameter :: acc_kind = merge(c_long_double, c_double,           &
       c_long_double > 0)

  real(acc_kind)               :: s

  s = 0.0_acc_kind
  do i = 1, n
     s = s + real(x(i), kind=acc_kind)
  end do

  value = real(s, kind=c_double)
end subroutine exprsum


      
