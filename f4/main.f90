program main
    use box
    implicit none
    interface
        subroutine write_matrix (a)
            real*4, dimension(:,:) :: a
        end subroutine write_matrix
        subroutine write_arr_int (a)
            integer*4, dimension(:) :: a
        end subroutine write_arr_int
        subroutine read_file(file_name, arr)
            character(LEN=*), intent(in) :: file_name
            real*4, dimension(:), allocatable, intent(out) :: arr
        end subroutine read_file
    end interface
    ! Grid
    integer*4, parameter :: ngrid = 256, ngx=128, ngy=64, ngz=256
    real*4, parameter :: side = 0.38, sx=0.3, sy=0.4, sz=0.2
    real*4, parameter, dimension(3) :: center = (/0.001,-0.001,0./)    
    real*4, dimension(:,:), allocatable :: grid
    
    ! Result
    real*4, dimension(:,:), allocatable :: boxed_extra
    
    ! Particles
    integer*4 :: npart, nextra, ninten, nexten, nvel
    real*4, dimension(:), allocatable :: arr_1D
    real*4, dimension(:,:), allocatable :: pos, vel
    real*4, dimension(:), allocatable :: mass, hsml, dens, pot, ene
    real*4, dimension(:,:), allocatable :: extra!, exten, inten
    
    ! Dummy
    integer*4 :: i
    
    PRINT '("Creando grilla...")'
    PRINT '(3X,"Celdas por lado: | ",3(I3," | "))', ngx, ngy, ngz
    PRINT '(3X,"Total de celdas: ",I9)', ngx * ngy * ngz
    PRINT '(3X,"Centro         : | ",3(F8.5," | "))', center
    PRINT '(3X,"Ancho de lados : | ",3(F8.5," | "))', sx, sy, sz
    PRINT '(2X,"Listo.")'
    call make_grid(ngx, ngy, ngz, sx, sy, sz, center, grid)
    PRINT '(1X,"Guardando grilla...")'
    OPEN(10, file="grilla", status="replace", form="unformatted", access="stream")
    WRITE(10) grid
    CLOSE(10)
    deallocate(grid)
    PRINT '(2X,"Listo.")'
    PRINT*, ""
    
    ! Leyendo posiciones
    PRINT '("Leyendo archivos...")'
    PRINT '(3X,"Leyendo posiciones...")'
    call read_file("pos_d", arr_1D)
    npart = SIZE(arr_1D) / 3
    allocate(pos(npart, 3))
    pos = RESHAPE(arr_1D, (/npart, 3/))
    deallocate(arr_1D)
    
    ! Leyendo el resto
    PRINT '(3X,"Leyendo masas...")'
    call read_file("mass_d", mass)
    PRINT '(3X,"Leyendo densidades...")'
    call read_file("dens_d", dens)
    PRINT '(3X,"Leyendo hsmls...")'
    call read_file("hsml_d", hsml)
    PRINT '(3X,"Leyendo potenciales...")'
    call read_file("pot_d", pot)
    PRINT '(3X,"Leyendo energías internas...")'
    call read_file("ene_d", ene)
    PRINT '(3X,"Leyendo velocidades...")'
    call read_file("vel_d", arr_1D)
    nvel = SIZE(arr_1D) / npart
    PRINT '(4X"Cantidad de velocidades: ",I1)', nvel
    allocate(vel(npart,nvel))
    vel = RESHAPE(arr_1D, (/npart, nvel/))
    PRINT '(2X,"Listo.")'
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
    PRINT '("Cantidad de valores totales   : ",I1)', nextra
    PRINT '("Cantidad de valores intensivos: ",I1)', ninten
    PRINT '("Cantidad de valores extensivos: ",I1)', nexten
    PRINT*, ""
    
    PRINT '("Calculando...")'
    call grid_extra(pos, mass, dens, hsml, extra, ninten,&
                    & ngx, ngy, ngz, sx, sy, sz, center,&
                    & boxed_extra)
    
    PRINT '("Escribiendo archivo de resultados ''valores'' con los valores grillados...")'
    OPEN(40, file="valores", status="replace", form="unformatted",access="stream")
    do i = 1, nextra
        WRITE(40) boxed_extra(i,:)
    end do
    CLOSE(40)
    PRINT '(2X,"Listo.")'
    PRINT '("Listo.")'
end program main

subroutine read_file(file_name, arr)
    character(LEN=*), intent(in) :: file_name
    real*4, dimension(:), allocatable, intent(out) :: arr
    logical :: existe
    integer*4 :: file_size, n_bytes
    real*4 :: this_byte
    INQUIRE(file=file_name, size=file_size, exist=existe)
    n_bytes = file_size / 4
    if (.not. existe) then
        PRINT*, "El archivo", file_name, "no existe."
        call EXIT(1)
    end if
    allocate(arr(n_bytes))
    OPEN(1, file=file_name, status="old", form="unformatted", access="stream")
    !! reading byte positions: 1 - 4, 4+1 - 4*2, ...
    do i = 1, n_bytes
        READ(1, pos=(i*4 - 3)) this_byte
        arr(i) = this_byte
    enddo   
    CLOSE(1)
end subroutine read_file

subroutine write_arr_int(a)
   integer*4, dimension(:) :: a
   WRITE(*,*)
   
   do i = LBOUND(a,1), UBOUND(a,1)
      WRITE(*,*) a(i)
   end do
end subroutine write_arr_int

subroutine write_matrix(a)
   real*4, dimension(:,:) :: a
   WRITE(*,*)
   
   do i = LBOUND(a,1), ubound(a,1)
      WRITE(*,*) (a(i,j), j = LBOUND(a,2), UBOUND(a,2))
   end do
end subroutine write_matrix
