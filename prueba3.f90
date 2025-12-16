module parameters
    implicit none
    integer, parameter :: N = 800, M = 20
    real(8), parameter :: dt = 0.000005d0, hx = 0.10d0, hy = 0.10d0
    real(8), parameter :: f = 3d0, q = 0.0002d0, e = 0.01d0
    real(8), parameter :: Du = 1d0, Dv = 1d0
    integer, parameter :: Tmax = 30000d0 * 1000d0
    real(8), parameter :: vl_initial = 0.98d0, vl_max =1.5d0, vl_step = 0.02d0
end module parameters

program CI

    use parameters
    implicit none

    ! Variables
    integer :: I, J, T, file_count, vl_step_count
    integer, parameter :: max_files = 900  ! Máximo de archivos por valor de vl
    double precision :: vl
    double precision, dimension(N, M) :: u, v, uf, vf
    double precision, dimension(N, M) :: u_initial, v_initial  ! Matrices para guardar los valores iniciales
    double precision, dimension(M) :: vel
    character(4) :: num
    character(20) :: folder_name

    ! Inicialización de variables
    vl = vl_initial
    vl_step_count = 0

    ! Leer el archivo CI.dat una vez y almacenar los valores iniciales en u_initial y v_initial
    call load_initial_conditions(u_initial, v_initial)

    ! Bucle principal para reiniciar el programa tras cada carpeta
    do while (vl <= vl_max)
        ! Crear carpeta para este valor de vl
        write(folder_name, '(A, F4.2)') 'v-', vl
        call system("mkdir "//trim(folder_name))

        ! Inicializar las condiciones con el nuevo valor de vl (usar las matrices iniciales)
        call initialize_from_constants(u, v, u_initial, v_initial, vel, vl)

        ! Bucle principal de simulación para este valor de vl
        file_count = 0
        do T = 1, Tmax
            call bz_reaction(u, v, uf, vf, vel)
            u = uf
            v = vf

            ! Aplicar condiciones de no flujo
            call no_flux_boundary(u, v)

            ! Aplicar condiciones de periodicidad
            call periodic_boundary(u, v)

            ! Guardar datos cada 5000 iteraciones
            if (mod(T, 3000) == 0) then
                write(num, '(I4.4)') T / 3000
                call save_data_to_file(N, M, u, v, trim(folder_name)//'/datosVel'//trim(num))

                ! Incrementar el contador de archivos
                file_count = file_count + 1

                ! Si se han generado 450 archivos, salir de este ciclo y aumentar vl
                if (file_count >= max_files) exit
            end if
        end do

        ! Incrementar el valor de vl para la siguiente iteración
        vl = vl + vl_step
    end do

    print *, "Simulación completada."
end program CI

! Subrutina para cargar las condiciones iniciales de CI.dat solo una vez
subroutine load_initial_conditions(u_initial, v_initial)
    use parameters
    implicit none
    double precision, dimension(N, M), intent(out) :: u_initial, v_initial
    integer :: II, JJ, III, JJJ, ios
    open(10, file='CI.dat', status='Old')

    ! Leer condiciones iniciales desde archivo y almacenarlas en las matrices u_initial y v_initial
    do II = 1, N
        do JJ = 1, M
            read(10, *, iostat=ios) III, JJJ, u_initial(III, JJJ), v_initial(III, JJJ)
            if (ios /= 0) exit  ! Si llega al final del archivo, salir del bucle
        end do
    end do

    print *, 'Condiciones iniciales cargadas: ', u_initial(10,10), v_initial(10,10)
    close(10)
end subroutine load_initial_conditions

! Subrutina para inicializar desde las matrices constantes
subroutine initialize_from_constants(u, v, u_initial, v_initial, vel, vl)
    use parameters
    implicit none
    double precision, dimension(N, M), intent(out) :: u, v
    double precision, dimension(N, M), intent(in) :: u_initial, v_initial
    double precision, dimension(M), intent(out) :: vel
    double precision, intent(in) :: vl  ! Pasar vl como argumento
    integer :: JJ
	real(8) :: h, y, v_max

    ! Asignar las condiciones iniciales desde las matrices constantes
    u = u_initial
    v = v_initial
    ! Definir la altura del canal y la velocidad máxima en el centro
    v_max = vl    ! Velocidad máxima en el centro
    ! Inicializar velocidades
	vel(1) =0
	vel(M-1) =0
    do JJ = 2, M-1  ! Coordenada y centrada en el canal
        vel(JJ) = v_max*JJ*(21-JJ) 
    end do

end subroutine initialize_from_constants

! Subrutina para guardar los datos en un archivo
subroutine save_data_to_file(n, m, u, v, filename)
    implicit none
    integer, intent(in) :: n, m
    double precision, dimension(n, m), intent(in) :: u, v
    character(*), intent(in) :: filename
    integer :: I, J
    logical :: valido

    ! Abrir archivo
    open(unit=40, file=trim(filename)//".dat", status='unknown')

    ! Escribir los datos de u y v para cada posición (i,j)
    do J = 1, m
        write(40, *) ' '
        do I = 1, n
            write(40, *) I, J, u(I, J), v(I, J)
        end do
    end do

    ! Cerrar archivo
    close(40)
end subroutine save_data_to_file

! Subrutina que realiza el cálculo de la reacción BZ
subroutine bz_reaction(u, v, uf, vf, vel)
    use parameters
    implicit none
    double precision, dimension(N, M), intent(in) :: u, v
    double precision, dimension(N, M), intent(out) :: uf, vf
    double precision, dimension(M), intent(in) :: vel
    integer :: I, J
    double precision :: Div, Dxu, Dyu, Dxv, Dyv, dxxu, dyyu, dxxv, dyyv

    ! Calcular la evolución de la reacción
    do I = 2, N-1
        do J = 2, M-1
            ! Primera ecuación
            Div = (u(I, J) - q) / (u(I, J) + q)
            Dxu = (u(I+1, J) + u(I-1, J) - 2*u(I, J)) / (hx**2)
            Dyu = (u(I, J+1) + u(I, J-1) - 2*u(I, J)) / (hy**2)
            dxxu = (u(I+1, J) - u(I-1, J)) / (2*hx)
            dyyu = (u(I, J+1) - u(I, J-1)) / (2*hy)
            uf(I, J) = u(I, J) + dt * (Du * (Dxu + Dyu) + (1/e) * (u(I, J) - u(I, J)**2 - f*v(I, J)*Div) - vel(J) * (dxxu))

            ! Segunda ecuación
            Dxv = (v(I+1, J) + v(I-1, J) - 2*v(I, J)) / (hx**2)
            Dyv = (v(I, J+1) + v(I, J-1) - 2*v(I, J)) / (hy**2)
            dxxv = (v(I+1, J) - v(I-1, J)) / (2*hx)
            dyyv = (v(I, J+1) - v(I, J-1)) / (2*hy)
            vf(I, J) = v(I, J) + dt * (Dv * (Dxv + Dyv) + u(I, J) - v(I, J) - vel(J) * (dxxv))
        end do
    end do
end subroutine bz_reaction

! Subrutina para aplicar las condiciones de no flujo
subroutine no_flux_boundary(u, v)
    use parameters
    implicit none
    double precision, dimension(N, M), intent(inout) :: u, v
    integer :: I

    ! No flujo en los bordes
    do I = 2, N-1
        u(I, 1) = u(I, 2)
        u(I, M) = u(I, M-1)
        v(I, 1) = v(I, 2)
        v(I, M) = v(I, M-1)
    end do
end subroutine no_flux_boundary

! Subrutina para aplicar las condiciones de periodicidad
subroutine periodic_boundary(u, v)
    use parameters
    implicit none
    double precision, dimension(N, M), intent(inout) :: u, v
    integer :: J

    ! Periodicidad en los extremos
    do J = 1, M
        u(1, J) = u(N-1, J)
        u(N, J) = u(2, J)
        v(1, J) = v(N-1, J)
        v(N, J) = v(2, J)
    end do
end subroutine periodic_boundary

