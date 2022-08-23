program main
    use omp_lib
    use box
    implicit none
    interface
        subroutine write_matrix (a)
            real, dimension(:,:) :: a
        end subroutine write_matrix
        subroutine write_arr_int (a)
            integer*4, dimension(:) :: a
        end subroutine write_arr_int
        subroutine read_file(file_name, nbytes, arr)
            character(LEN=*), intent(in) :: file_name
            integer*4, intent(in) :: nbytes
            real, dimension(:), allocatable, intent(out) :: arr
        end subroutine read_file
    end interface

    ! Grid
    integer*4, parameter :: ngrid = 128, ngx=256, ngy=256, ngz=256
    real, parameter :: side = 0.2165, sx=0.5, sy=0.5, sz=0.5
    real, parameter, dimension(3) :: center = (/0.,0.,0./)
    real, dimension(:,:), allocatable :: grid

    ! Result
    real, dimension(:,:), allocatable :: boxed_extra

    ! Particles
    integer*4 :: npart, nextra, ninten, nexten, nvel
    real, dimension(:), allocatable :: arr_1D_p, arr_1D_v
    real, dimension(:,:), allocatable :: pos, vel
    real, dimension(:), allocatable :: mass, hsml, dens, pot, ene
    real, dimension(:,:), allocatable :: extra!, exten, inten
    integer*4, parameter :: nbytes_in = 4 ! Los datos están en real*8 o real*4 ?

    ! Dummy
    integer*4 :: i
    real*8 :: time, time2

    ! Parallel
    integer*4 :: nthreads = 4
    nthreads = MAX(MIN(nthreads, OMP_GET_MAX_THREADS()-2),1)
    call OMP_SET_NUM_THREADS(nthreads)

    PRINT '("Trabajando con ",I2, " threads.")', OMP_GET_MAX_THREADS()
    PRINT '("Las operaciones (y archivos de salida) serán en ",I1, " bytes.")', SIZEOF(side)
    PRINT '("Leyendo archivos se leerán en ",I1, " bytes.")', nbytes_in
    !$OMP PARALLEL IF((nthreads<=4).AND.(1<nthreads)) DEFAULT(PRIVATE) &
    !$OMP SHARED(mass,dens,hsml,pot,ene,nvel,vel,pos,npart,arr_1D_p,arr_1D_v)
    !$OMP SECTIONS 
    !$OMP SECTION
        PRINT '(3X,"Leyendo posiciones...")'
        call read_file("pos_d", nbytes_in, arr_1D_p)
    !$OMP SECTION
        PRINT '(3X,"Leyendo masas...")'
        call read_file("mass_d", nbytes_in, mass)
    !$OMP SECTION
        PRINT '(3X,"Leyendo densidades...")'
        call read_file("dens_d", nbytes_in, dens)
    !$OMP SECTION
        PRINT '(3X,"Leyendo hsmls...")'
        call read_file("hsml_d", nbytes_in, hsml)
    !$OMP SECTION
        PRINT '(3X,"Leyendo potenciales...")'
        call read_file("pot_d", nbytes_in, pot)
    !$OMP SECTION
        PRINT '(3X,"Leyendo energías internas...")'
        call read_file("ene_d", nbytes_in, ene)
    !$OMP SECTION
        PRINT '(3X,"Leyendo velocidades...")'
        call read_file("vel_d", nbytes_in, arr_1D_v)
    !$OMP END SECTIONS
    !$OMP END PARALLEL
    !! Parámetros...
    PRINT '(2X,"Obteniendo parámetros...")'
    npart = SIZE(arr_1D_p) / 3
    allocate(pos(npart, 3))
    pos = RESHAPE(arr_1D_p, (/npart, 3/))
    deallocate(arr_1D_p)
    nvel = SIZE(arr_1D_v) / npart
    allocate(vel(npart,nvel))
    vel = RESHAPE(arr_1D_v, (/npart, nvel/))
    deallocate(arr_1D_v)
    PRINT '(1X,"Listo.")'
    PRINT*, ""

    ! Juntamos
    !! Intensivos: dens, vel(x3), pot
    !!Extensivos: masa, ene
    nextra = 4 + nvel
    allocate(extra(npart, nextra))
    extra(:,1) = dens
    do i = 2, 1 + nvel
        extra(:,i) = vel(:, i - 1)
    end do
    extra(:, 2 + nvel) = pot
    extra(:, 3 + nvel) = mass
    extra(:, 4 + nvel) = ene

    nexten = 2
    ninten = nextra - nexten
    PRINT '("Cantidad de partículas: ",I9)', npart
    PRINT '("Cantidad de velocidades: ",I1)', nvel
    PRINT '("Cantidad de valores totales   : ",I1)', nextra
    PRINT '("Cantidad de valores intensivos: ",I1)', ninten
    PRINT '("Cantidad de valores extensivos: ",I1)', nexten
    PRINT*, ""

    time = OMP_GET_WTIME()
    PRINT '("Generando grilla...")'
    PRINT '(3X,"Celdas por lado: | ",3(I3," | "))', ngx, ngy, ngz
    PRINT '(3X,"Total de celdas: ",I9)', ngx * ngy * ngz
    PRINT '(3X,"Centro         : | ",3(F8.5," | "))', center
    PRINT '(3X,"Ancho de lados : | ",3(F8.5," | "))', sx, sy, sz
    call make_grid(ngx, ngy, ngz, sx, sy, sz, center, grid)
    PRINT '(2X,"Guardando grilla...")'
    OPEN(10, file="grilla", status="replace", form="unformatted", access="stream")
    WRITE(10) grid
    CLOSE(10)
    deallocate(grid)
    PRINT '(2X,"Tiempo de generación y guardado: ",F8.4, " [s].")', OMP_GET_WTIME() - time
    PRINT '(1X,"Listo.")'
    PRINT*, ""

    PRINT '("Calculando...")'
    call grid_extra(pos, mass, dens, hsml, extra, ninten,&
                    & ngx, ngy, ngz, sx, sy, sz, center,&
                    & boxed_extra)
   

    time2 = OMP_GET_WTIME()
    PRINT '("Escribiendo archivo de resultados ''valores'' con los valores grillados...")'
    OPEN(40, file="valores", status="replace", form="unformatted",access="stream")
    do i = 1, nextra
        WRITE(40) boxed_extra(i,:)
    end do
    CLOSE(40)
    PRINT '("Tiempo de guardado de valores: ",F8.4, " [s].")', OMP_GET_WTIME() - time2
    PRINT '(2X,"Listo.")'
    PRINT '("Tiempo TOTAL: ",F8.4, " [s].")', OMP_GET_WTIME() - time
    PRINT '("Listo.")'
end program main

subroutine read_file(file_name, nbytes, arr)
    use omp_lib
    character(LEN=*), intent(in) :: file_name
    integer*4, intent(in) :: nbytes
    real, dimension(:), allocatable, intent(out) :: arr
    logical :: existe
    integer*4 :: file_size, n_bytes, i, thread
    real*4 :: this_byte_4
    real*8 :: this_byte_8

    INQUIRE(file=file_name, size=file_size, exist=existe)
    n_bytes = file_size / nbytes
    if (.not. existe) then
        PRINT*, "El archivo", file_name, "no existe."
        call EXIT(1)
    end if
    allocate(arr(n_bytes))
    thread = OMP_GET_THREAD_NUM()
    OPEN(thread, file=file_name, status="old", form="unformatted", access="stream")
    !! reading byte positions: 1 - nbytes, nbytes+1 - nbytes*2, ...
    if (nbytes == 8) then
        do i = 1, n_bytes
            READ(thread, pos=(i*nbytes - nbytes + 1)) this_byte_8
            arr(i) = REAL(this_byte_8, SIZEOF(arr(i)))
        enddo
    else
        do i = 1, n_bytes
            READ(thread, pos=(i*nbytes - nbytes + 1)) this_byte_4
            arr(i) = REAL(this_byte_4, SIZEOF(arr(i)))
        enddo
    end if
    CLOSE(thread)
end subroutine read_file

subroutine write_arr_int(a)
   integer*4, dimension(:) :: a
   WRITE(*,*)
   do i = LBOUND(a,1), UBOUND(a,1)
      WRITE(*,*) a(i)
   end do
end subroutine write_arr_int

subroutine write_matrix(a)
   real, dimension(:,:) :: a
   WRITE(*,*)
   do i = LBOUND(a,1), ubound(a,1)
      WRITE(*,*) (a(i,j), j = LBOUND(a,2), UBOUND(a,2))
   end do
end subroutine write_matrix
