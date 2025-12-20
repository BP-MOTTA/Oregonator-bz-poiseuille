module parameters
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: N = 500, M = 21
    real(dp), parameter :: hx = 0.10_dp, hy = 0.10_dp
    real(dp), parameter :: dt = 5.0e-6_dp
    real(dp), parameter :: f = 3.0_dp, q = 2.0e-4_dp, e = 1.0e-2_dp
    real(dp), parameter :: Du = 1.0_dp, Dv = 1.0_dp
    integer, parameter :: Tmax = 30000*1000

    real(dp), parameter :: vl_initial = 0.00_dp
    real(dp), parameter :: vl_max     = 1.50_dp
    real(dp), parameter :: vl_step    = 0.5_dp

    integer, parameter :: Ndiag = 5000
    integer, parameter :: Nsnap = 1500
end module parameters


program CI
    use parameters
    use io_hdf5
    implicit none

    integer :: T, file_count
    integer, parameter :: max_files = 1800          ! solo para .dat
    real(dp) :: vl
    real(dp), dimension(N, M) :: u, v, uf, vf
    real(dp), dimension(N, M) :: u_initial, v_initial
    real(dp), dimension(M)    :: vel

    character(4)   :: num
    character(32)  :: folder_name
    character(256) :: h5name
    integer :: vl_int, vl_whole, vl_frac

    integer :: Nx_phys
    real(dp), allocatable :: us(:,:), vs(:,:)

    ! --- setup general ---
    vl = vl_initial
    Nx_phys = N - 2
    allocate(us(Nx_phys, M), vs(Nx_phys, M))

    call load_initial_conditions(u_initial, v_initial)

    do while (file_count <= 9000)

        ! --- nombre de carpeta vX_YY ---
        vl_int   = nint(vl * 100.0_dp)
        vl_whole = vl_int / 100
        vl_frac  = abs(mod(vl_int, 100))
        write(folder_name, '(A,I0,"_",I2.2)') 'v', vl_whole, vl_frac

        call system('mkdir -p '//trim(folder_name))
        print *, "Run folder => ", trim(folder_name)

        call initialize_from_constants(u, v, u_initial, v_initial, vel, vl)

        call iniciar_diag_xcm(trim(folder_name)//'/xcm_x.dat')

        h5name = trim(folder_name)//'/run.h5'
        print *, "HDF5 output => ", trim(h5name)
        call h5_init_run(h5name, Nx_phys, M, hx, hy, dt, vl, N, M)

        file_count = 0

        do T = 1, Tmax

            call bz_reaction(u, v, uf, vf, vel)
            u = uf
            v = vf

            call no_flux_boundary(u, v)
            call periodic_boundary(u, v)

            ! Diagnóstico xcm
            if (mod(T, Ndiag) == 0) then
                call diagnostico_xcm_x(T, u, trim(folder_name)//'/xcm_x.dat')
            end if

            ! Snapshot HDF5
            if (mod(T, Nsnap) == 0) then
                us(:,:) = u(2:N-1, 1:M)
                vs(:,:) = v(2:N-1, 1:M)
                call h5_write_snapshot(dt*real(T,dp), us, vs)
                file_count=file_count +1
            end if

            ! Salida .dat (opcional / pesado)
            if (mod(T, 300000) == 0) then
                write(num, '(I4.4)') file_count
                call save_data_to_file(N, M, u, v, trim(folder_name)//'/datosVel'//trim(num))
                if (file_count >= max_files) exit
            end if

        end do

        call h5_close_run()

        vl = vl + vl_step
    end do

    deallocate(us, vs)
    print *, 'Simulación completada.'
end program CI


subroutine load_initial_conditions(u_initial, v_initial)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(out) :: u_initial, v_initial
    integer :: i_read, j_read, ios
    real(dp) :: u_val, v_val

    ! importante: evitar basura si el archivo no llena todo
    u_initial = 0.0_dp
    v_initial = 0.0_dp

    open(10, file='CI.dat', status='old', action='read')
    do
        read(10, *, iostat=ios) i_read, j_read, u_val, v_val
        if (ios /= 0) exit
        if (i_read >= 1 .and. i_read <= N .and. j_read >= 1 .and. j_read <= M) then
            u_initial(i_read, j_read) = u_val
            v_initial(i_read, j_read) = v_val
        end if
    end do
    close(10)
    print *, 'Condiciones iniciales cargadas desde CI.dat'
end subroutine load_initial_conditions


subroutine initialize_from_constants(u, v, u_initial, v_initial, vel, vl)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(out) :: u, v
    real(dp), dimension(N, M), intent(in)  :: u_initial, v_initial
    real(dp), dimension(M),    intent(out) :: vel
    real(dp), intent(in) :: vl
    integer :: j
    real(dp) :: y

    u = u_initial
    v = v_initial

    do j = 1, M
        y = real(j-1,dp) / real(M-1,dp)
        vel(j) = vl * 4.0_dp * y * (1.0_dp - y)   ! vmax = vl
    end do
end subroutine initialize_from_constants


subroutine save_data_to_file(n_in, m_in, u, v, filename)
    use parameters, only: dp
    implicit none
    integer, intent(in) :: n_in, m_in
    real(dp), dimension(n_in, m_in), intent(in) :: u, v
    character(*), intent(in) :: filename
    integer :: i, j

    open(unit=40, file=trim(filename)//'.dat', status='replace', action='write')
    do j = 1, m_in
        write(40,*) ' '
        do i = 1, n_in
            write(40,*) i, j, u(i,j), v(i,j)
        end do
    end do
    close(40)
end subroutine save_data_to_file


subroutine bz_reaction(u, v, uf, vf, vel)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(in)  :: u, v
    real(dp), dimension(N, M), intent(out) :: uf, vf
    real(dp), dimension(M),    intent(in)  :: vel
    integer :: i, j
    real(dp) :: div, Dxu, Dyu, Dxv, Dyv
    real(dp) :: ux_adv, vx_adv

    !$omp parallel do collapse(2) default(shared) private(i,j,div,Dxu,Dyu,Dxv,Dyv,ux_adv,vx_adv) schedule(static)
    do i = 2, N-1
        do j = 2, M-1

            div = (u(i,j) - q) / (u(i,j) + q)

            Dxu = (u(i+1,j) + u(i-1,j) - 2.0_dp*u(i,j)) / (hx*hx)
            Dyu = (u(i,j+1) + u(i,j-1) - 2.0_dp*u(i,j)) / (hy*hy)

            if (vel(j) >= 0.0_dp) then
                ux_adv = (u(i,j) - u(i-1,j)) / hx
            else
                ux_adv = (u(i+1,j) - u(i,j)) / hx
            end if

            uf(i,j) = u(i,j) + dt * ( Du*(Dxu + Dyu) &
                 + (1.0_dp/e) * (u(i,j) - u(i,j)**2 - f*v(i,j)*div) &
                 - vel(j) * ux_adv )

            Dxv = (v(i+1,j) + v(i-1,j) - 2.0_dp*v(i,j)) / (hx*hx)
            Dyv = (v(i,j+1) + v(i,j-1) - 2.0_dp*v(i,j)) / (hy*hy)

            if (vel(j) >= 0.0_dp) then
                vx_adv = (v(i,j) - v(i-1,j)) / hx
            else
                vx_adv = (v(i+1,j) - v(i,j)) / hx
            end if

            vf(i,j) = v(i,j) + dt * ( Dv*(Dxv + Dyv) + u(i,j) - v(i,j) &
                 - vel(j) * vx_adv )
        end do
    end do
    !$omp end parallel do

    ! bordes: se ajustan con BCs luego
    uf(1,:) = u(1,:);  uf(N,:) = u(N,:)
    uf(:,1) = u(:,1);  uf(:,M) = u(:,M)

    vf(1,:) = v(1,:);  vf(N,:) = v(N,:)
    vf(:,1) = v(:,1);  vf(:,M) = v(:,M)
end subroutine bz_reaction


subroutine no_flux_boundary(u, v)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(inout) :: u, v
    integer :: i

    do i = 2, N-1
        u(i,1) = u(i,2);     u(i,M) = u(i,M-1)
        v(i,1) = v(i,2);     v(i,M) = v(i,M-1)
    end do
end subroutine no_flux_boundary


subroutine periodic_boundary(u, v)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(inout) :: u, v
    integer :: j

    do j = 1, M
        u(1,j) = u(N-1,j)
        u(N,j) = u(2,j)
        v(1,j) = v(N-1,j)
        v(N,j) = v(2,j)
    end do
end subroutine periodic_boundary


subroutine iniciar_diag_xcm(filename)
    use parameters
    implicit none
    character(len=*), intent(in) :: filename
    integer :: unit

    unit = 60
    open(unit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(A)') '# t  x_cm  theta_cm(rad)  umax  flag'
    close(unit)
end subroutine iniciar_diag_xcm


subroutine diagnostico_xcm_x(t, u, filename)
    use parameters
    implicit none
    integer, intent(in) :: t
    real(dp), dimension(N,M), intent(in) :: u
    character(len=*), intent(in) :: filename

    integer :: i, j, unit, flag
    integer :: Nx_phys
    real(dp) :: mass, umax, t_phys
    real(dp) :: theta, csum, ssum, uij
    real(dp) :: theta_cm, xcm
    real(dp) :: pi, Lx

    pi = acos(-1.0_dp)

    Nx_phys = N - 2
    Lx = real(Nx_phys,dp) * hx

    mass = 0.0_dp
    umax = 0.0_dp
    csum = 0.0_dp
    ssum = 0.0_dp

    do i = 2, N-1
        theta = 2.0_dp*pi * real(i-2,dp) / real(Nx_phys,dp)
        do j = 1, M
            uij  = u(i,j)
            mass = mass + uij
            csum = csum + uij * cos(theta)
            ssum = ssum + uij * sin(theta)
            if (uij > umax) umax = uij
        end do
    end do

    if (umax > 1.0e-3_dp) then
        flag = 1
    else
        flag = 0
    end if

    if (mass > 0.0_dp) then
        theta_cm = atan2(ssum, csum)
        if (theta_cm < 0.0_dp) theta_cm = theta_cm + 2.0_dp*pi
        xcm = (theta_cm / (2.0_dp*pi)) * Lx
    else
        theta_cm = -1.0_dp
        xcm      = -1.0_dp
    end if

    t_phys = dt * real(t,dp)

    unit = 60
    open(unit=unit, file=trim(filename), status='unknown', position='append', action='write')
    write(unit,'(F15.8,1X,F15.8,1X,F15.8,1X,E15.8,1X,I2)') t_phys, xcm, theta_cm, umax, flag
    close(unit)
end subroutine diagnostico_xcm_x
