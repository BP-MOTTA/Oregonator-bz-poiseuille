module parameters
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: N = 400, M = 21 !seleccionar adecuadamente segun la malla, en este caso es de 400X21 para la simetria del Pouseille
    real(dp), parameter :: hx = 0.10_dp, hy = 0.10_dp
    real(dp), parameter :: dt = 5.0e-6_dp
    real(dp), parameter :: f = 3.0_dp, q = 2.0e-4_dp, e = 1.0e-2_dp !parametros de estabilidad del BZ (Leer en literatura)
    real(dp), parameter :: Du = 1.0_dp, Dv = 1.0_dp !Difusion de u y v
    integer, parameter :: Tmax = 30000*1000      ! 3.0e7 pasos
    real(dp), parameter :: vl_initial = 0.00_dp
    real(dp), parameter :: vl_max     = 1.50_dp
    real(dp), parameter :: vl_step    = 0.02_dp
    integer, parameter :: Ndiag = 5000 !lectura de datos de velocidad
end module parameters


program CI
    use parameters
    implicit none

    ! Variables principales
    integer :: I, J, T, file_count
    integer, parameter :: max_files = 900  ! máx. archivos por vl
    real(dp) :: vl !valor de velocidad del flujo promedio
    real(dp), dimension(N, M) :: u, v, uf, vf !variables de concentracion quimica
    real(dp), dimension(N, M) :: u_initial, v_initial !variables de condiciones iniciales
    real(dp), dimension(M)    :: vel        !velocidad del flujo laminar
    character(4)  :: num
    character(32) :: folder_name
    integer :: vl_int, vl_whole, vl_frac

    vl = vl_initial
    call load_initial_conditions(u_initial, v_initial) !lectura de las CI


    do while (vl <= vl_max + 1.0e-10_dp) !cambio de los valores promedios de velocida
        vl_int   = nint(vl * 100.0_dp)          !hundredths as integer
        vl_whole = vl_int / 100
        vl_frac  = abs(mod(vl_int, 100))
        write(folder_name, '(A,I0,"_",I2.2)') 'v', vl_whole, vl_frac
        call system('mkdir -p '//trim(folder_name)) 
        call initialize_from_constants(u, v, u_initial, v_initial, vel, vl) 
        print*, folder_name
        call iniciar_diag_xcm(trim(folder_name)//'/xcm_x.dat')
        ! Bucle de tiempo
        file_count = 0
        do T = 1, Tmax

            call bz_reaction(u, v, uf, vf, vel)
            u = uf
            v = vf

            call no_flux_boundary(u, v)
            call periodic_boundary(u, v)
            
            ! Guardar datos cada 3000 pasos
            if (mod(T, 3000) == 0) then
                write(num, '(I4.4)') T / 300000
                call save_data_to_file(N, M, u, v, &
                     trim(folder_name)//'/datosVel'//trim(num))

                file_count = file_count + 1
                if (file_count >= max_files) exit
            end if
            !lectura de velocidad
            if (mod(T, Ndiag) == 0) call diagnostico_xcm_x(T, u, trim(folder_name)//'/xcm_x.dat')

        end do

        ! Siguiente valor de velocidad máxima
        vl = vl + vl_step
    end do

    print *, 'Simulación completada.'
end program CI

!cargado de las CI desde los archivos generados por CI
subroutine load_initial_conditions(u_initial, v_initial)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(out) :: u_initial, v_initial
    integer :: i_read, j_read
    real(dp) :: u_val, v_val
    integer :: ios
    open(10, file='CI.dat', status='old', action='read')
    ! CI.dat/ I,J,u(i,j), v(i,j)
    do
        read(10, *, iostat=ios) i_read, j_read, u_val, v_val !lectura directa
        if (ios /= 0) exit !sistema antierrores de lectura
        if (i_read >= 1 .and. i_read <= N .and. j_read >= 1 .and. j_read <= M) then
            u_initial(i_read, j_read) = u_val
            v_initial(i_read, j_read) = v_val
        end if
    end do

    close(10)
    print *, 'Condiciones iniciales cargadas desde CI.dat'
end subroutine load_initial_conditions


!============================================================
!  Inicializar campos + perfil de velocidad de Poiseuille
!============================================================
subroutine initialize_from_constants(u, v, u_initial, v_initial, vel, vl)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(out) :: u, v
    real(dp), dimension(N, M), intent(in)  :: u_initial, v_initial
    real(dp), dimension(M),    intent(out) :: vel
    real(dp), intent(in) :: vl
    integer :: j
    real(dp) :: y

    ! Copiar CI
    u = u_initial
    v = v_initial

    ! Perfil de Poiseuille:
    ! y* = (j-1)/(M-1) en [0,1], v(y) = vl * 4 y (1-y)
    do j = 1, M
        y = real(j-1,dp) / real(M-1,dp)
        vel(j) = vl * 4.0_dp * y * (1.0_dp - y)
    end do
end subroutine initialize_from_constants


!============================================================
!  Guardar datos en archivo
!============================================================
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


!============================================================
!  Paso de reacción–difusión–advección BZ + Poiseuille
!============================================================
subroutine bz_reaction(u, v, uf, vf, vel)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(in)  :: u, v
    real(dp), dimension(N, M), intent(out) :: uf, vf
    real(dp), dimension(M),    intent(in)  :: vel
    integer :: i, j
    real(dp) :: div, Dxu, Dyu, Dxv, Dyv
    real(dp) :: ux_adv, vx_adv

    ! Esquema explícito, difusión central + advección upwind
    !$omp parallel do collapse(2) private(div, Dxu, Dyu, Dxv, Dyv, ux_adv, vx_adv, i, j)
    do i = 2, N-1
        do j = 2, M-1

            ! --- Reacción + difusión en u ---
            div = (u(i,j) - q) / (u(i,j) + q)
            Dxu = (u(i+1,j) + u(i-1,j) - 2.0_dp*u(i,j)) / (hx*hx)
            Dyu = (u(i,j+1) + u(i,j-1) - 2.0_dp*u(i,j)) / (hy*hy)

            ! Término advectivo en x: vel(j) * du/dx, upwind
            if (vel(j) >= 0.0_dp) then
                ux_adv = (u(i,j) - u(i-1,j)) / hx
            else
                ux_adv = (u(i+1,j) - u(i,j)) / hx
            end if

            uf(i,j) = u(i,j) + dt * ( Du*(Dxu + Dyu) &
                 + (1.0_dp/e) * (u(i,j) - u(i,j)**2 - f*v(i,j)*div) &
                 - vel(j) * ux_adv )

            ! --- Reacción + difusión en v ---
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

    ! Copiamos bordes tal cual; se corrigen con las BCs
    uf(1,:) = u(1,:)
    uf(N,:) = u(N,:)
    uf(:,1) = u(:,1)
    uf(:,M) = u(:,M)

    vf(1,:) = v(1,:)
    vf(N,:) = v(N,:)
    vf(:,1) = v(:,1)
    vf(:,M) = v(:,M)
end subroutine bz_reaction


!============================================================
!  Condiciones de no flujo en y
!============================================================
subroutine no_flux_boundary(u, v)
    use parameters
    implicit none
    real(dp), dimension(N, M), intent(inout) :: u, v
    integer :: i

    do i = 2, N-1
        u(i,1) = u(i,2)
        u(i,M) = u(i,M-1)
        v(i,1) = v(i,2)
        v(i,M) = v(i,M-1)
    end do
end subroutine no_flux_boundary


!============================================================
!  Condiciones periódicas en x
!============================================================
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
    ! t  x_cm  theta_cm  umax  flag
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

    ! Dominio físico en x: i = 2 .. N-1  (porque 1 y N son ghost)
    Nx_phys = N - 2
    Lx = real(Nx_phys,dp) * hx

    mass = 0.0_dp
    umax = 0.0_dp
    csum = 0.0_dp
    ssum = 0.0_dp

    do i = 2, N-1
        ! i=2 -> theta=0, i=N-1 -> theta=2pi*(Nx_phys-1)/Nx_phys
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
        theta_cm = atan2(ssum, csum)      ! [-pi, pi]
        if (theta_cm < 0.0_dp) theta_cm = theta_cm + 2.0_dp*pi  ! [0,2pi)
        xcm = (theta_cm / (2.0_dp*pi)) * Lx                      ! [0,Lx)
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