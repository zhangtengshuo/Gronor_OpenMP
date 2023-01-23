subroutine normalize(a,n)
implicit none

integer,intent(in)         :: n
integer                    :: i
real(kind=8),intent(inout) :: a(n)
real(kind=8)               :: norm

norm = 0.0d0
do i = 1, n
  norm = norm + a(i) * a(i)
end do
norm = 1.0d0 / sqrt(norm)
do i = 1, n
  a(i) = a(i) * norm
end do

return
end subroutine normalize
