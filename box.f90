module box
    use kern
    implicit none

    contains

            integer*4 function ijk2boxid(i, j, k, ngrid) result (id)
                implicit none
                integer*4, intent(in) :: i, j, k, ngrid

                if ((i < 0) .or. (i >= ngrid)) then
                    id = -1
                else if ((j < 0) .or. (j >= ngrid)) then
                    id = -1
                else if ((k < 0) .or. (k >= ngrid)) then
                    id = -1
                else
                    id = i + j * ngrid + k * ngrid * ngrid
                end if
            end function ijk2boxid !0 -> ngrid³ -1
            
            integer*4 function boxid2k(id, ngrid) result (k)
                implicit none
                integer*4, intent(in) :: id, ngrid

                k = id / (ngrid * ngrid)
            end function boxid2k !0 -> ngrid³ -1
            
            integer*4 function boxid2j(id, ngrid) result (j)
                implicit none
                integer*4, intent(in) :: id, ngrid

                j = MOD(id / ngrid, ngrid)
            end function boxid2j
            
            integer*4 function boxid2i(id, ngrid) result (i)
                implicit none
                integer*4, intent(in) :: id, ngrid

                i = MOD(id, ngrid)
            end function boxid2i
            
            subroutine boxid2ijk (id, ngrid, i, j, k)
                implicit none
                integer*4, intent(in) :: id, ngrid
                integer*4, intent(out) :: i, j, k

                i = boxid2i(id, ngrid)
                j = boxid2j(id, ngrid)
                k = boxid2k(id, ngrid)
            end subroutine boxid2ijk

            subroutine grid_pos (pos, npart, ngrid, side, center, box_id)
                implicit none
                real*8, dimension(npart,3), intent(in) :: pos
                integer*4, intent(in) :: ngrid, npart
                real*8, intent(in) :: side
                real*8, dimension(3), intent(in) :: center
                integer*4, dimension(npart), intent(inout) :: box_id
                real*8, dimension(npart,3) :: centred, pos_norm
                integer*4, dimension(npart,3) :: box_ijk
                logical, dimension(npart) :: in_grid
                real*8 :: h, side_2
                integer*4 :: i, j

                h = side / (1. * ngrid)
                side_2 = side * 0.5
                in_grid = .true.
                do j = 1, 3
                    do i = 1, npart
                        centred(i,j) = pos(i,j) - center(j) ! centramos las partículas
                        ! llevamos al cuadrante positivo y normalizamos en h
                        pos_norm(i,j) = (centred(i,j) + (0.5 * side)) / (1. * h) 
                        in_grid(i) = in_grid(i) .and. (ABS(centred(i,j)) < side_2) ! cuáles están dentro de la grilla?
                        box_ijk(i, j) = FLOOR(pos_norm(i,j)) ! cuál es el box (índice izquierdo)
                    end do
                end do
                do i = 1, npart ! cuáles están dentro de la grilla?
                    if (in_grid(i)) then
                        box_id(i) = ijk2boxid(box_ijk(i,1), box_ijk(i,2), box_ijk(i,3), ngrid)
                    else 
                        box_id(i) = -1
                    end if
                end do

            end subroutine grid_pos

            subroutine grid_extra (pos, hsml, mass, extra, support, npart, ngrid, side, center, boxed_data)
                implicit none
                real*8, dimension(npart,3), intent(in) :: pos
                real*8, dimension(npart), intent(in) :: hsml, mass
                real*8, dimension(:, :), intent(in) :: extra ! (npart, n_extra)
                real*8, intent(in) :: support, side
                integer*4, intent(in) :: npart, ngrid                
                real*8, dimension(3), intent(in) :: center
                real*8, dimension(:,:), allocatable, intent(out) :: boxed_data

                integer*4, dimension(npart) :: box_id
                real*8, dimension(npart,3) :: centred, pos_norm
                integer*4, dimension(npart,3) :: box_ijk
                logical, dimension(npart) :: in_grid

                real*8 :: h, side_2, d2prev, d2next, d, W
                integer*4 :: i, j, p, n_useful, ngrid3, n_extra
                integer*4 :: bi, bj, bk, bin

                real*8, allocatable, dimension(:) :: hsml_norm_u, mass_u, vol_norm_u
                real*8, allocatable, dimension(:,:) :: pos_norm_u
                integer*4, allocatable, dimension(:,:) :: box_ijk_u, first, last
                real*8, allocatable, dimension(:) :: Sj, Ik
                real*8, allocatable, dimension(:,:) :: Aj ! este va a ser (n_extra, npart)
                real*8, allocatable, dimension(:,:) :: Ak ! este va a ser (n_extra, ng³)

                if (SIZE(extra, 1) .ne. npart) then
                    WRITE(*,*) "Bad 'extra' input shape:"
                    WRITE(*,*) " ", SIZE(extra, 1), " was given, and ", npart, " was expected"
                    WRITE(*,*) "Exiting."
                    call exit(1)
                end if

                h = side / (1. * ngrid) ! ancho del grid
                side_2 = side * 0.5  ! mitad de ancho del grid
                in_grid = .true.
                do j = 1, 3 ! cols
                    do i = 1, npart ! rows
                        centred(i,j) = pos(i,j) - center(j) ! centramos las partículas
                        pos_norm(i,j) = (centred(i,j) + (0.5 * side)) / h ! llevamos al cuadrante positivo y normalizamos en h
                        in_grid(i) = in_grid(i) .and. (ABS(centred(i,j)) < side_2) ! cuáles están dentro de la grilla?
                        box_ijk(i, j) = FLOOR(pos_norm(i,j)) ! cuál es el box (índice izquierdo)
                    end do
                end do
                do i = 1, npart 
                    if (in_grid(i)) then ! cuáles están dentro de la grilla?
                        box_id(i) = ijk2boxid(box_ijk(i, 1), box_ijk(i, 2), box_ijk(i, 3), ngrid) !va de 0 a ngrid-1
                    else 
                        box_id(i) = -1
                    end if
                end do
                WRITE(*,*) "Índice de boxes de partículas hallados."
                if ((support <= 0) .or. ANY(hsml < 0)) then ! si no hay más nada, ya fue
                    return
                end if
                n_useful = COUNT(in_grid) ! esta cantidad usaremos
                ! alocatamos
                allocate(hsml_norm_u(n_useful)) 
                allocate(pos_norm_u(n_useful,3), box_ijk_u(n_useful,3))
                ! ponemos en estos array los datos correspondientes
                hsml_norm_u = PACK(hsml * support / h, in_grid) ! multiplicamos y normalizamos en h
                mass_u = PACK(mass, in_grid)
                vol_norm_u = 4. * PI * hsml_norm_u * hsml_norm_u * hsml_norm_u / 3.
                do j = 1, 3
                    pos_norm_u(:,j) = PACK(pos_norm(:,j), in_grid)
                    box_ijk_u(:,j) = PACK(box_ijk(:,j), in_grid)
                end do

                ! ahora vemos cuáles box tocan
                allocate(first(n_useful,3), last(n_useful,3))
                do j = 1, 3
                    do i = 1, n_useful
                        d2prev = MOD(pos_norm_u(i,j), 1.)
                        d2next = 1. - d2prev
                        first(i,j) = MAX(box_ijk_u(i,j) - FLOOR(ABS(hsml_norm_u(i) - d2prev), kind=4), 0)
                        last(i,j) = MIN(box_ijk_u(i,j) + FLOOR(ABS(hsml_norm_u(i) - d2next), kind=4), ngrid-1)
                    end do
                end do
                WRITE(*,*) "Índice de boxes vecinos de partículas hallados."
                
                WRITE(*,*) "Calculando pesos Sj..."
                ! Ahora calculamos Sj
                allocate(Sj(n_useful))
                Sj = 0.
                do i = 1, n_useful
                    do bk = first(i,3), last(i,3)
                        do bj = first(i,2), last(i,2)
                            do bi = first(i,1), last(i,1)
                                d = SUM((pos_norm_u(i,:) - (/bi, bj, bk/) - 0.5)**2)
                                Sj(i) = Sj(i) + spline_d(d, hsml_norm_u(i))
                            end do
                        end do
                    end do
                end do
                WRITE(*,*) "Sj calculados."
                WRITE(*,*) "Calculando pesos Ak..."
                ! Ahora calculamos Ak
                n_extra = SIZE(extra, 2)
                ngrid3 = ngrid * ngrid * ngrid
                allocate(Ak(n_extra,ngrid3), Ik(ngrid3))
                allocate(Aj(n_extra,n_useful))
                do j = 1, n_extra
                    Aj(j,:) = PACK(extra(:,j), in_grid)
                end do
                Aj = 1.
                do i = 1, n_useful
                    if (Sj(i) > 0.) then
                        do bk = first(i,3), last(i,3)
                            do bj = first(i,2), last(i,2)
                                do bi = first(i,1), last(i,1)
                                    bin = ijk2boxid(bi, bj, bk, ngrid) + 1
                                    d = SUM((pos_norm_u(i,:) - (/bi, bj, bk/) - 0.5)**2)
                                    W = spline_d(d, hsml_norm_u(i))
                                    do p = 1, n_extra
                                        Ak(p,bin) = Ak(p,bin) + (vol_norm_u(i) / Sj(i) * W * Aj(p,i))
                                    end do
                                    Ik(bin) = Ik(bin) + (vol_norm_u(i) * W)
                                end do
                            end do
                        end do
                    else
                        d = SUM((pos_norm_u(i,:) - (/box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3)/) - 0.5)**2)
                        W = spline_d(d, hsml_norm_u(i))
                        bin = ijk2boxid(box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3), ngrid) + 1
                        do p = 1, n_extra
                            Ak(p,bin) = Ak(p,bin) + (vol_norm_u(i) * Aj(p,i))
                        end do
                        Ik(bin) = Ik(bin) + (vol_norm_u(i) * W)
                    end if
                end do
                deallocate(hsml_norm_u, box_ijk_u)
                deallocate(pos_norm_u, mass_u, vol_norm_u)
                deallocate(Sj, Aj)
                deallocate(first, last)
                WRITE(*,*) "Ak calculados."
                do i = 1, ngrid3
                    if (Ik(i) > 0) then
                        Ak(:,i) = Ak(:,i) / Ik(i)
                    end if
                end do
                allocate(boxed_data(n_extra,ngrid3))
                boxed_data = Ak
                deallocate(Ak, Ik)
            end subroutine grid_extra
    
end module box
