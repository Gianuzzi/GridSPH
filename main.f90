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
    ! Particles
    !! Parameters
!     integer*4, parameter :: npart_des = 5000000
!     integer*4, parameter :: nside = FLOOR(npart_des**(1/3.))
!     integer*4, parameter :: npart = nside * nside * nside
    !! Arrays
!     real*8, dimension(npart) :: mass, hsml
!     real*8, dimension(npart,2) :: extra
!     real*8, dimension(:,:), allocatable :: pos
    
    ! Grid
    integer*4, parameter :: ngrid = 256
    real*8, parameter :: side = 0.38
    real*8, parameter, dimension(3) :: center = (/0.,0.,0./)    
    real*8, dimension(:,:), allocatable :: grid
    
    ! Result
    integer*4, dimension(:), allocatable :: box_id
    real*8, dimension(:,:), allocatable :: boxed_extra
    
    ! V2. Particles
    integer*4 :: npart, nextra
    real*8, dimension(:), allocatable :: arr_1D
    real*8, dimension(:,:), allocatable :: pos
    real*8, dimension(:), allocatable :: mass, hsml, dens
    real*8, dimension(:,:), allocatable :: extra
    
    ! Dummy
    integer*4 :: i
    
    PRINT*, "Creando grilla con ngrid...", ngrid
    PRINT*, "  ngrid³:", ngrid**3
    PRINT*, "  centro:", center
    call make_grid(ngrid, side, center, grid)
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
    PRINT*, "Cantidad de partículas:", npart
    allocate(pos(npart,3))
    pos = RESHAPE(arr_1D, (/npart,3/))
    deallocate(arr_1D)
    
    ! Leyendo masa, densidad y hsml
    PRINT*, "Leyendo masa..."
    call read_file("mass_d", mass)
    PRINT*, "Leyendo densidad..."
    call read_file("dens_d", dens)
    PRINT*, "Leyendo hsml..."
    call read_file("hsml_d", hsml)
    ! Leyendo "extra"
!     nextra = 1
!     allocate(extra(npart,nextra))
!     extra(:,1) = mass
    PRINT*, "Leyendo extra..."
    call read_file("vel_d", arr_1D)
    nextra = SIZE(arr_1D) / npart
    PRINT*, "Cantidad de datos extra:", nextra
    allocate(extra(npart,nextra))
    extra = RESHAPE(arr_1D, (/npart,nextra/))
    
    PRINT*, "Calculando..."
    allocate(box_id(npart))
    call grid_extra(pos, mass, dens, hsml, extra, ngrid, side, center, box_id, boxed_extra)
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
