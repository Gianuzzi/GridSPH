module box
    use kern
    implicit none

    contains

            integer*4 function ijk2boxid(i, j, k, ngrid) result (id)
                !(i, j, k, NGRID) -> box_id
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
                !(box_id, NGRID) -> k
                implicit none
                integer*4, intent(in) :: id, ngrid

                k = id / (ngrid * ngrid)
            end function boxid2k !0 -> ngrid³ -1
            
            integer*4 function boxid2j(id, ngrid) result (j)
                !(box_id, NGRID) -> j
                implicit none
                integer*4, intent(in) :: id, ngrid

                j = MOD(id / ngrid, ngrid)
            end function boxid2j
            
            integer*4 function boxid2i(id, ngrid) result (i)
                !(box_id, NGRID) -> i
                implicit none
                integer*4, intent(in) :: id, ngrid

                i = MOD(id, ngrid)
            end function boxid2i
            
            subroutine boxid2ijk (id, ngrid, i, j, k)
                !from (box_id, NGRID) -> set: (i, j, k)
                implicit none
                integer*4, intent(in) :: id, ngrid
                integer*4, intent(out) :: i, j, k

                i = boxid2i(id, ngrid)
                j = boxid2j(id, ngrid)
                k = boxid2k(id, ngrid)
            end subroutine boxid2ijk
            
            subroutine make_grid (ngrid, side, center, grid)
                implicit none
                integer*4, intent(in) :: ngrid
                real*8, intent(in) :: side
                real*8, dimension(3), intent(in) :: center
                real*8, dimension(:,:), allocatable, intent(out) :: grid
                integer*4 :: i, j, k, m
                real*8 :: ratio 
                
                ratio = side / (ngrid * 1.)                
                allocate(grid(ngrid * ngrid * ngrid, 3))
                do k = 1, ngrid
                    do j = 1, ngrid
                        do i = 1, ngrid
                            m = ijk2boxid(i-1, j-1, k-1, ngrid) + 1
                            grid(m,3) = k
                            grid(m,2) = j
                            grid(m,1) = i
                        end do
                    end do
                end do
                grid = grid * ratio
                do m = 1, ngrid * ngrid * ngrid
                    grid(m,:) = grid(m,:) - 0.5 * (side + ratio) * (/1, 1, 1/) + center
                end do
            end subroutine make_grid

            subroutine grid_extra (pos, mass, dens, hsml, extra, ngrid, side, center, box_id, boxed_extra)
                implicit none
                real*8, dimension(:), intent(in) :: mass
                real*8, dimension(SIZE(mass)), intent(in) :: dens, hsml
                real*8, dimension(SIZE(mass),3), intent(in) :: pos
                real*8, dimension(:, :), intent(in) :: extra ! (npart, n_extra)
                real*8, intent(in) :: side
                integer*4, intent(in) :: ngrid                
                real*8, dimension(3), intent(in) :: center
                real*8, dimension(:,:), allocatable, intent(out) :: boxed_extra
                integer*4, dimension(SIZE(mass)), intent(out) :: box_id
                
                real*8, dimension(SIZE(mass),3) :: centred, pos_norm
                integer*4, dimension(SIZE(mass),3) :: box_ijk
                logical, dimension(SIZE(mass)) :: in_grid

                real*8 :: h, side_2, d2prev, d2next, d, W
                integer*4 :: i, j, p, n_useful, ngrid3, n_extra
                integer*4 :: bi, bj, bk, bin

                real*8, allocatable, dimension(:) :: mass_u, dens_norm_u, hsml_norm_u, vol_norm_u
                real*8, allocatable, dimension(:,:) :: pos_norm_u
                integer*4, allocatable, dimension(:,:) :: box_ijk_u, first, last
                real*8, allocatable, dimension(:) :: Sj, Ik
                real*8, allocatable, dimension(:,:) :: Aj ! este va a ser (n_extra, npart)
                real*8, allocatable, dimension(:,:) :: Ak ! este va a ser (n_extra, ng³)

                ! Checkeamos que "extra" esté bien
                if (SIZE(extra, 1) .ne. SIZE(mass)) then
                    WRITE(*,*) "Bad 'extra' input shape:"
                    WRITE(*,*) " ", SIZE(extra, 1), " was given, and ", SIZE(mass), " was expected"
                    WRITE(*,*) "Exiting."
                    call EXIT(1)
                end if

                ! Calculamos parámetros básicos
                h = side / (1. * ngrid) ! ancho de cada box del grid
                side_2 = side * 0.5  ! mitad de ancho del grid
                ngrid3 = ngrid * ngrid * ngrid ! esta es la cantidad total de boxes
                n_extra = SIZE(extra, 2) ! esta es la cantidad de propiedades a calcular

                ! Normalizamos las posiciones y vemos cuáles partículas caen dentro del grid
                in_grid = .true.
                do j = 1, 3 ! columnas
                    do i = 1, SIZE(mass) ! filas
                        centred(i,j) = pos(i,j) - center(j) ! centramos las partículas
                        pos_norm(i,j) = (centred(i,j) + (0.5 * side)) / h ! llevamos al "cuadrante positivo" y normalizamos en h
                        in_grid(i) = in_grid(i) .and. (ABS(centred(i,j)) < side_2) ! in_grid son las que caen dentro del grid
                        box_ijk(i, j) = FLOOR(pos_norm(i,j)) ! calculamos el box(i,j,k) de cada partícula. De 0 a ngrid - 1
                    end do
                end do
                ! Pasamos de (i,j,k) a box_id
                do i = 1, SIZE(mass) 
                    if (in_grid(i)) then ! las que están dentro de la grilla...
                        box_id(i) = ijk2boxid(box_ijk(i, 1), box_ijk(i, 2), box_ijk(i, 3), ngrid) !va de 0 a ngrid-1
                    else 
                        box_id(i) = -1 !-1 tienen las que caen fuera de la grilla
                    end if
                end do
                WRITE(*,*) "Índice de boxes de partículas hallados."
                
                ! Cantidad de partículas que caen dentro de la grilla
                n_useful = COUNT(in_grid) ! esta es la cantidad de partículas ("útiles") que usaremos
                PRINT*, "Cantidad de partículas dentro del grid:", n_useful
                PRINT*, "útiles/total = ", n_useful / (1. * SIZE(mass))
                
                ! Alocatamos 
                allocate(mass_u(n_useful), dens_norm_u(n_useful), hsml_norm_u(n_useful), vol_norm_u(n_useful)) 
                allocate(pos_norm_u(n_useful,3), box_ijk_u(n_useful,3))
                allocate(first(n_useful,3), last(n_useful,3))
                allocate(Sj(n_useful))
                allocate(Aj(n_extra,n_useful))
                allocate(Ak(n_extra,ngrid3), Ik(ngrid3))

                ! Ponemos en estos array los datos correspondientes. Todos normalizados en h=lado_de_box
                hsml_norm_u = PACK(hsml, in_grid) / h ! hsml normalizado en h
                mass_u = PACK(mass, in_grid) ! masa
                dens_norm_u = PACK(dens, in_grid) * (h * h * h) ! densidad normalizada
                vol_norm_u = mass_u / dens_norm_u
!                 vol_norm_u = 4. * PI * hsml_norm_u * hsml_norm_u * hsml_norm_u / 3. ! volumen normalizado
!                 vol_norm_u = 1. ! ESTO PARA EXTRÍNSECOS
                do j = 1, 3
                    pos_norm_u(:,j) = PACK(pos_norm(:,j), in_grid) ! posición normalizada
                    box_ijk_u(:,j) = PACK(box_ijk(:,j), in_grid) ! (i,j,k) del box. De 0 a ngrid - 1
                end do
                do j = 1, n_extra
                    Aj(j,:) = PACK(extra(:,j), in_grid) ! cantidad "extra" (sin normalizar!!)
                end do

                ! Calculamos a cuáles box toca cada partícula, según su hsml y posición
                do j = 1, 3
                    do i = 1, n_useful
                        d2prev = MOD(pos_norm_u(i,j), 1.) ! que tan metida está en su box
                        d2next = 1. - d2prev ! cuanto le falta para el box siguiente
                        first(i,j) = box_ijk_u(i,j) - FLOOR(ABS(1 + hsml_norm_u(i) - d2prev), kind=4) ! primero
                        last(i,j) = box_ijk_u(i,j) + FLOOR(ABS(1 + hsml_norm_u(i) - d2next), kind=4) ! último
                    end do
                end do
                WRITE(*,*) "Índice de boxes vecinos de partículas hallados."
                
                WRITE(*,*) "Calculando pesos Sj..."
                ! Calculamos los pesos Sj de cada partícula
                Sj = 0.
                do i = 1, n_useful
                    do bk = first(i,3), last(i,3) ! del primero al último
                        do bj = first(i,2), last(i,2) ! del primero al último
                            do bi = first(i,1), last(i,1) ! del primero al último
                                d = SQRT(SUM((pos_norm_u(i,:) - (/bi, bj, bk/) - 0.5)**2)) ! desde su posición al centro del otro box
                                Sj(i) = Sj(i) + spline_d(d, hsml_norm_u(i)) ! sumamos contribución
                            end do
                        end do
                    end do
                end do
                WRITE(*,*) "Sj calculados."
                OPEN(70, file="Sj", status="replace", form="unformatted",access="stream")
                WRITE(70) Sj
                CLOSE(70)
                
                ! Arreglamos First y Last
                do j = 1, 3
                    do i = 1, n_useful
                        first(i,j) = MAX(MIN(first(i,j), ngrid - 1), 0) ! primero
                        last(i,j) = MIN(MAX(last(i,j), 0), ngrid - 1) ! último
                    end do
                end do
                
                WRITE(*,*) "Calculando pesos Ak..."
                ! Calculamos los pesos Ak de cada box. Acá hay que tener en cuenta el valor de Aj ("extra") de cada partícula
                ! También calculamos las constantes de normalización Ik
                do i = 1, n_useful
                    if (Sj(i) > 0.) then ! en este caso debemos recorrer los boxes amigos
                        do bk = first(i,3), last(i,3) ! del primero al último
                            do bj = first(i,2), last(i,2) ! del primero al último
                                do bi = first(i,1), last(i,1) ! del primero al último
                                    bin = ijk2boxid(bi, bj, bk, ngrid) + 1
                                    d = SQRT(SUM((pos_norm_u(i,:) - (/bi, bj, bk/) - 0.5)**2))
                                    W = spline_d(d, hsml_norm_u(i))
                                    do p = 1, n_extra
                                        Ak(p,bin) = Ak(p,bin) + (vol_norm_u(i) / Sj(i) * W * Aj(p,i))
                                    end do
                                    Ik(bin) = Ik(bin) + (vol_norm_u(i) * W) ! cte de normalización 
                                end do
                            end do
                        end do
                    else ! en este caso, solo consideramos el box más cercano (el suyo)
                        d = SQRT(SUM((pos_norm_u(i,:) - (/box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3)/) - 0.5)**2))
                        W = spline_d(d, hsml_norm_u(i))
                        bin = ijk2boxid(box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3), ngrid) + 1
                        do p = 1, n_extra
                            Ak(p,bin) = Ak(p,bin) + (vol_norm_u(i) * Aj(p,i))
                        end do
                        Ik(bin) = Ik(bin) + (vol_norm_u(i) * W) ! cte de normalización
                    end if
                end do
                ! Dealocatamos
                deallocate(hsml_norm_u, box_ijk_u)
                deallocate(pos_norm_u, mass_u, vol_norm_u)
                deallocate(Sj, Aj)
                deallocate(first, last)
                WRITE(*,*) "Ak calculados."
                ! Normalizamos según Ik (SOLO SIRVE PARA INTRÍNSECOS)
                do i = 1, ngrid3
                    if (Ik(i) > 0) then
                        Ak(:,i) = Ak(:,i) / Ik(i)
                    end if
                end do
                ! Obtemos resultado final
                allocate(boxed_extra(n_extra,ngrid3))
                boxed_extra = Ak
                deallocate(Ak, Ik)
            end subroutine grid_extra
    
end module box
