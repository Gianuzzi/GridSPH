program main
    use box
    implicit none
    interface
        subroutine write_matrix (a)
            real*8, dimension(:,:) :: a
        end subroutine write_matrix
        subroutine write_arr_int (a)
            integer*4, dimension(:) :: a
        end subroutine write_arr_int
        subroutine read_file(file_name, arr)
            character(LEN=*), intent(in) :: file_name
            real*8, dimension(:), allocatable, intent(out) :: arr
        end subroutine read_file
    end interface
    ! Grid
    integer*4, parameter :: ngrid = 256, ngx=128, ngy=128, ngz=128
    real*8, parameter :: side = 0.38, sx=0.2, sy=0.2, sz=0.2
    real*8, parameter, dimension(3) :: center = (/0.,0.,0./)    
    real*8, dimension(:,:), allocatable :: grid
    
    ! Result
    integer*4, dimension(:), allocatable :: box_id
    real*8, dimension(:,:), allocatable :: boxed_extra
    
    ! Particles
    integer*4 :: npart, nextra, ninten, nexten, nvel
    real*8, dimension(:), allocatable :: arr_1D
    real*8, dimension(:,:), allocatable :: pos, vel
    real*8, dimension(:), allocatable :: mass, hsml, dens, pot, ene
    real*8, dimension(:,:), allocatable :: extra!, exten, inten
    
    ! Dummy
    integer*4 :: i
    
    PRINT*, "Creando grilla con:"
    PRINT*, "  Celdas por lado:", ngx, ngy, ngz
    PRINT*, "  Total de celdas:", ngx * ngy * ngz
    PRINT*, "  Centro:", center
    PRINT*, "  Ancho de lados:", sx, sy, sz
    call make_grid(ngx, ngy, ngz, sx, sy, sz, center, grid)
    PRINT*, "Guardando grilla..."
    OPEN(10, file="grilla", status="replace", form="unformatted", access="stream")
    WRITE(10) grid
    CLOSE(10)
    deallocate(grid)
    PRINT*, "Listo."
    
    ! Leyendo posiciones
    PRINT*, "Leyendo posiciones..."
    call read_file("pos_d", arr_1D)
    npart = SIZE(arr_1D) / 3
    PRINT*, " Cantidad de partículas:", npart
    allocate(pos(npart, 3))
    pos = RESHAPE(arr_1D, (/npart, 3/))
    deallocate(arr_1D)
    
    ! Leyendo masa, densidad y hsml
    PRINT*, "Leyendo mas..."
    call read_file("mass_d", mass)
    PRINT*, "Leyendo densidad..."
    call read_file("dens_d", dens)
    PRINT*, "Leyendo hsml..."
    call read_file("hsml_d", hsml)
    PRINT*, "Leyendo potencial..."
    call read_file("pot_d", pot)
    PRINT*, "Leyendo energia interna..."
    call read_file("ene_d", ene)
    
    PRINT*, "Leyendo velocidad..."
    call read_file("vel_d", arr_1D)
    nvel = SIZE(arr_1D) / npart
    PRINT*, " Cantidad de velocidades:", nvel
    allocate(vel(npart,nvel))
    vel = RESHAPE(arr_1D, (/npart, nvel/))
    
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
    PRINT*, "Cantidad de valores totales   :", nextra
    PRINT*, "Cantidad de valores intensivos:", ninten
    PRINT*, "Cantidad de valores extensivos:", nexten
    
    PRINT*, "Calculando..."
    allocate(box_id(npart))
    call grid_extra(pos, mass, dens, hsml, extra, ninten,&
                    & ngx, ngy, ngz, sx, sy, sz, center,&
                    & box_id, boxed_extra)
    PRINT*, "Guardando resultados..."
    PRINT*, " box_id."
    OPEN(30, file="box_id", status="replace", form="unformatted",access="stream")
    WRITE(30) box_id
    CLOSE(30)
    PRINT*, " valores."
    OPEN(40, file="valores", status="replace", form="unformatted",access="stream")
    do i = 1, nextra
        WRITE(40) boxed_extra(i,:)
    end do
    CLOSE(40)
    PRINT*, "Listo."
end program main

subroutine read_file(file_name, arr)
    character(LEN=*), intent(in) :: file_name
    real*8, dimension(:), allocatable, intent(out) :: arr
    logical :: existe
    integer*4 :: file_size, n_bytes
    real*8 :: this_byte
    INQUIRE(file=file_name, size=file_size, exist=existe)
    n_bytes = file_size / 8
    if (.not. existe) then
        PRINT*, "El archivo", file_name, "no existe."
        call EXIT(1)
    end if
    allocate(arr(n_bytes))
    OPEN(1, file=file_name, status="old", form="unformatted", access="stream")
    !! reading byte positions: 1 - 8, 8+1 - 8*2, ...
    do i = 1, n_bytes
        READ(1, pos=(i*8 - 7)) this_byte
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
   real*8, dimension(:,:) :: a
   WRITE(*,*)
   
   do i = LBOUND(a,1), ubound(a,1)
      WRITE(*,*) (a(i,j), j = LBOUND(a,2), UBOUND(a,2))
   end do
end subroutine write_matrix
