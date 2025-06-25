module parametros
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: N=400, M=21
    real(dp),parameter :: dt = 1.0d-6, hx = 0.1d0, hy = 0.1d0
    real(dp), parameter :: f = 3.0d0, q = 0.0002d0, e = 0.01d0
    real(dp), parameter :: Du = 1.0d0, Dv = 1.0d0
    integer, parameter :: Tmax = 30000 * 1000
end module parametros


module variables
  use parametros
  implicit none
  real(dp), dimension(N, M) :: u, v, uf, vf
end module variables

program CI_modular
  use parametros
  use variables
  implicit none
  integer :: i, j, t
  real(dp) :: itst4
  character(4) :: num

  call inicializar()
  call aplicar_pulso()

  do t = 1, Tmax
    call paso_bz()
    call condiciones_frontera()
    if (t == 300000) call eliminar_pulso()
    itst4 = mod(t, 5000)
    if (itst4 == 0) then
      write(num, '(I4.4)') t / 5000
      call guardar_texto("CI_salida/" // num // ".dat")
    end if
  end do
end program CI_modular

subroutine inicializar()
  use variables
  integer :: i, j
  do i = 1, N
    do j = 1, M
      u(i,j) = 0.00039984_dp
      v(i,j) = 0.00039984_dp
    end do
  end do
end subroutine

subroutine aplicar_pulso()
  use variables
  integer :: i, j
  do i = N/2 - 5, N/2 + 5
    do j = 1, M
      u(i,j) = 50 * 0.00039984_dp
    end do
  end do
end subroutine

subroutine paso_bz()
  use parametros
  use variables
  integer :: i, j
  real(dp) :: div, dxu, dyu, dxv, dyv
  do i = 2, N - 1
    do j = 2, M - 1
      div = (u(i,j) - q) / (u(i,j) + q)
      dxu = (u(i+1,j) + u(i-1,j) - 2*u(i,j)) / (hx**2)
      dyu = (u(i,j+1) + u(i,j-1) - 2*u(i,j)) / (hy**2)
      uf(i,j) = u(i,j) + dt * (Du * (dxu + dyu) + (1/e) * (u(i,j) - u(i,j)**2 - f*v(i,j)*div))

      dxv = (v(i+1,j) + v(i-1,j) - 2*v(i,j)) / (hx**2)
      dyv = (v(i,j+1) + v(i,j-1) - 2*v(i,j)) / (hy**2)
      vf(i,j) = v(i,j) + dt * (Dv * (dxv + dyv) + u(i,j) - v(i,j))
    end do
  end do
  u = uf; v = vf
end subroutine

subroutine condiciones_frontera()
  use variables
  integer :: i, j
  do i = 2, N - 1
    u(i,1) = u(i,2); u(i,M) = u(i,M-1)
    v(i,1) = v(i,2); v(i,M) = v(i,M-1)
  end do
  do j = 1, M
    u(1,j) = u(N-1,j); u(N,j) = u(2,j)
    v(1,j) = v(N-1,j); v(N,j) = v(2,j)
  end do
end subroutine

subroutine eliminar_pulso()
  use variables
  integer :: i, j
  do i = 1, N/2
    do j = 1, M
      u(i,j) = 0.00039984_dp
      v(i,j) = 0.00039984_dp
    end do
  end do
end subroutine

subroutine guardar_texto(filename)
  use variables
  character(*) :: filename
  integer :: i, j, unit
  unit = 40
  open(unit=unit, file=filename, status='replace')
  do j = 1, M
    write(unit,*) ' '
    do i = 1, N
      write(unit,*) i, j, u(i,j), v(i,j)
    end do
  end do
  close(unit)
end subroutine

