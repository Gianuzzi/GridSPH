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
    end interface
    integer*4, parameter :: npart_des = 10
    integer*4, parameter :: ngrid = 256
    integer*4, parameter :: nside = FLOOR(npart_des**(1/3.))
    integer*4, parameter :: npart = nside * nside * nside
    real*8, parameter :: side = 1.
    real*8, parameter :: support = 2.
    real*8, parameter, dimension(3) :: center = (/0.,0.,0./)
    real*8, dimension(:,:), allocatable :: boxed_data
    real*8, parameter :: ratio = side / (nside * 1.)
    integer *4 :: i, j, k, m
    real*8, dimension(npart,3) :: pos
    real*8, dimension(npart) :: hsml, mass
    real*8, dimension(npart,2) :: extra
    print*,"Cantidad total de partículas:", npart
    do k = 1, nside
        do j = 1, nside
            do i = 1, nside
                m = (i - 1) + (j - 1) * nside + (k - 1) * nside * nside + 1
                pos(m,3) = (k - 1)
                pos(m,2) = (j - 1)
                pos(m,1) = (i - 1)
            end do
        end do
    end do
    pos = pos * ratio

    do m = 1, npart
        pos(m,:) = pos(m,:) - 0.5 * (side - ratio) * (/1, 1, 1/)
    end do
    call write_matrix(pos)
    mass = 1.
    hsml = 0.1 * ratio
    extra(:,1) = mass
    extra(:,2) = 2 * mass
    call grid_extra(pos, hsml, mass, extra, support, npart, ngrid, side, center, boxed_data)
    if (ANY(boxed_data < 0)) then
        WRITE(*,*) "Menor a 0"
    end if
    ! call write_matrix(TRANSPOSE(boxed_data))
end program main

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