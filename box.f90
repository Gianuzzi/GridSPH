module box
    use kern
    use omp_lib
    implicit none

    contains
            
        integer*4 function ijk2boxid(i, j, k, nx, ny, nz) result (id)
            !(i, j, k, nx, ny, nz) -> box_id
            implicit none
            integer*4, intent(in) :: i, j, k, nx, ny, nz

            if ((i < 0) .or. (i >= nx)) then
                id = -1
            else if ((j < 0) .or. (j >= ny)) then
                id = -1
            else if ((k < 0) .or. (k >= nz)) then
                id = -1
            else
                id = i + j * nx + k * nx * ny
            end if
        end function ijk2boxid !0 -> ngrid³ -1
        
        integer*4 function boxid2k(id, nx, ny) result (k)
            !(box_id, nx, ny) -> k
            implicit none
            integer*4, intent(in) :: id, nx, ny

            k = id / (nx * ny)
        end function boxid2k !0 -> nz - 1
        
        integer*4 function boxid2j(id, nx, ny) result (j)
            !(box_id, nx, ny) -> j
            implicit none
            integer*4, intent(in) :: id, nx, ny

            j = MOD(id / nx, ny)
        end function boxid2j !0 -> ny - 1
        
        integer*4 function boxid2i(id, nx) result (i)
            !(box_id, nx) -> i
            implicit none
            integer*4, intent(in) :: id, nx

            i = MOD(id, nx)
        end function boxid2i !0 -> nx - 1
        
        subroutine boxid2ijk (id, nx, ny, i, j, k)
            !from (box_id, nx, ny) -> set: (i, j, k)
            implicit none
            integer*4, intent(in) :: id, nx, ny
            integer*4, intent(out) :: i, j, k

            i = boxid2i(id, nx)
            j = boxid2j(id, nx, ny)
            k = boxid2k(id, nx, ny)
        end subroutine boxid2ijk !0,0,0 -> nx-1,ny-1,nz-1
        
        subroutine make_grid (nx, ny, nz, sx, sy, sz, center, grid)
            implicit none
            integer*4, intent(in) :: nx, ny, nz
            real, intent(in) :: sx, sy, sz
            real, dimension(3), intent(in) :: center
            real, dimension(:,:), allocatable, intent(out) :: grid
            integer*4 :: i, j, k, m
            real :: rx, ry , rz
            
            rx = sx / (nx * 1.)
            ry = sy / (ny * 1.)
            rz = sz / (nz * 1.)
            allocate(grid(nx * ny * nz, 3))
            !$OMP PARALLEL PRIVATE (k,j,i,m) FIRSTPRIVATE(nx,ny,nz,rx,ry,rz) SHARED (grid)
            !$OMP DO SCHEDULE (STATIC)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            m = ijk2boxid(i-1, j-1, k-1, nx, ny, nz) + 1
                            ! grid(m,1) = i * rx
                            ! grid(m,2) = j * ry
                            ! grid(m,3) = k * rz
                            grid(m,:) = (/i * rx, j * ry, k * rz/) - 0.5 * (/sx + rx, sy + ry, sz + rz/) + center
                        end do
                    end do
                end do
            !$OMP END DO
            !$OMP END PARALLEL
            ! do m = 1, nx * ny * nz
            !     grid(m,:) = grid(m,:) - 0.5 * (/sx + rx, sy + ry, sz + rz/) + center
            ! end do
        end subroutine make_grid
                
        subroutine grid_extra (pos, mass, dens, hsml, extra, ninten, ngx, ngy, ngz, sx, sy, sz, center, boxed_extra)
            implicit none
            real, dimension(:), intent(in) :: mass
            real, dimension(SIZE(mass)), intent(in) :: dens, hsml
            real, dimension(SIZE(mass),3), intent(in) :: pos
            real, dimension(:, :), intent(in) :: extra ! (npart, n_extra)
            real, intent(in) :: sx, sy, sz
            integer*4, intent(in) :: ngx, ngy, ngz, ninten
            real, dimension(3), intent(in) :: center
            real, dimension(:,:), allocatable, intent(out) :: boxed_extra
            
            real, dimension(SIZE(mass),3) :: centred, pos_norm
            integer*4, dimension(SIZE(mass),3) :: box_ijk
            integer*4, dimension(SIZE(mass)) :: box_id
            logical, dimension(SIZE(mass)) :: in_grid, touch_grid

            real, dimension(3) :: h, s, s2, h2
            real*8 :: start_time, time_1, time_2
            integer*4, dimension(3) :: ng
            real :: hx, hy, hz, sx2, sy2, sz2, d2prev, d2next, d, W, vol_k, dx2, dy2, dz2
            integer*4 :: i, j, n_useful, ngrid3, n_extra, nexten
            integer*4 :: bi, bj, bk, bin
            integer*4, allocatable, dimension(:) :: box_id_u
            real, allocatable, dimension(:) :: mass_u, dens_u, hsml_u, vol_u
            real, allocatable, dimension(:,:) :: pos_norm_u
            integer*4, allocatable, dimension(:,:) :: box_ijk_u, first, last
            real, allocatable, dimension(:) :: Sj, Ik
            real, allocatable, dimension(:,:) :: Aj,Aj_Sj ! este va a ser (n_extra, nuseful)
            real, allocatable, dimension(:,:) :: Ak ! este va a ser (n_extra, ng³)

            ! integer*4, allocatable, dimension(:) :: npart_per_cell, acum_part !(ng³)
            ! integer*4, allocatable, dimension(:) :: id_part_in_cell ! (nuseful)
            ! integer*4 :: prev
            
            
            ! Checkeamos que "extra" esté bien
            if (SIZE(extra, 1) .ne. SIZE(mass)) then
                PRINT*, "Bad 'extra' input shape:"
                PRINT*, " ", SIZE(extra, 1), " was given, and ", SIZE(mass), " was expected"
                PRINT*, "Exiting."
                call EXIT(1)
            end if
            if (SIZE(extra, 2) .lt. ninten) then
                PRINT*, "Bad 'extra' input shape:"
                PRINT*, " ", SIZE(extra, 2), " was given, and not lower than", ninten, " was expected"
                PRINT*, "Exiting."
                call EXIT(1)
            end if

            ! Calculamos parámetros básicos
            hx = sx / (1. * ngx) ! ancho de cada box del grid
            hy = sy / (1. * ngy) ! ancho de cada box del grid
            hz = sz / (1. * ngz) ! ancho de cada box del grid
            sx2 = sx * 0.5  ! mitad de ancho del grid
            sy2 = sy * 0.5  ! mitad de ancho del grid
            sz2 = sz * 0.5  ! mitad de ancho del grid
            ngrid3 = ngx * ngy * ngz ! esta es la cantidad total de boxes
            vol_k = hx * hy * hz ! volumen de las celdas
            s = (/sx, sy, sz/)
            s2 = (/sx2, sy2, sz2/)
            h = (/hx, hy, hz/)
            h2 = h * h
            ng = (/ngx, ngy, ngz/)
            n_extra = SIZE(extra, 2) ! esta es la cantidad de propiedades a calcular
            nexten = n_extra - ninten

            start_time = OMP_GET_WTIME()
            time_1 = OMP_GET_WTIME()
            
            PRINT*, ""
            PRINT*, "------ Comenzando subroutina 'grid_extra' ------"
            PRINT*, ""
            
            PRINT*, "Normalizando posiciones y detectando partículas interactuantes..."
            ! Normalizamos las posiciones y vemos cuáles partículas caen dentro del grid
            in_grid = .TRUE.
            touch_grid = .TRUE.
            do j = 1, 3 ! columnas
                do i = 1, SIZE(mass) ! filas
                    centred(i,j) = pos(i,j) - center(j) ! centramos las partículas
                    pos_norm(i,j) = (centred(i,j) + s2(j)) / h(j) ! llevamos al "cuadrante positivo" y normalizamos en h
                    in_grid(i) = in_grid(i) .and. (ABS(centred(i,j)) < s2(j)) ! in_grid son las que caen dentro del grid
                    touch_grid(i) = touch_grid(i) .and. ((ABS(centred(i,j)) - hsml(i)) < s2(j)) ! touch_grid son las que caen dentro del grid
                    box_ijk(i, j) = FLOOR(pos_norm(i,j)) ! calculamos el box(i,j,k) de cada partícula. De 0 a ng() - 1
                end do
            end do
            
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Calculando índice de box de partículas..."
            ! Pasamos de (i,j,k) a box_id
            do i = 1, SIZE(mass) 
                if (in_grid(i)) then ! las que están dentro de la grilla...
                    box_id(i) = ijk2boxid(box_ijk(i, 1), box_ijk(i, 2), box_ijk(i, 3), ngx, ngy, ngz) !va de 0 a ng-1
                else 
                    box_id(i) = -1 !-1 tienen las que caen fuera de la grilla
                end if
            end do
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Escribiendo archivo 'box_id' con índice de box de partículas..."
            OPEN(90, file="box_id", status="replace", form="unformatted", access="stream")
            WRITE(90) box_id
            CLOSE(90)
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()

            PRINT*, "Reduciendo set de datos..."
            ! Cantidad de partículas que caen dentro de la grilla
            ! !touch_grid = in_grid
            n_useful = COUNT(touch_grid) ! esta es la cantidad de partículas ("útiles") que usaremos
            PRINT '(3X,"Partículas dentro del grid           : ",I10)', COUNT(in_grid)
            PRINT '(3X,"Dentro/Total                         : ",F6.4)', COUNT(in_grid) / (1. * SIZE(mass))
            PRINT '(3X,"Partículas interactuantes con el grid: ",I10)', n_useful
            PRINT '(3X,"Dentro/Interactuantes                : ",F6.4)', COUNT(in_grid) / (1. * n_useful)
            PRINT '(3X,"Interactuantes/Total                 : ",F6.4)', n_useful / (1. * SIZE(mass))
            
            ! Alocatamos
            allocate(box_id_u(n_useful)) 
            allocate(mass_u(n_useful), dens_u(n_useful), hsml_u(n_useful), vol_u(n_useful)) 
            allocate(pos_norm_u(n_useful,3), box_ijk_u(n_useful,3))
            allocate(first(n_useful,3), last(n_useful,3))
            allocate(Sj(n_useful))
            allocate(Aj(n_extra,n_useful))
            allocate(Ak(n_extra,ngrid3), Ik(ngrid3))
            allocate(Aj_Sj(n_extra,n_useful))
            ! allocate(npart_per_cell(ngrid3), acum_part(ngrid3))
            ! allocate(id_part_in_cell(n_useful))
            
            ! Ponemos en estos array los datos correspondientes.
            box_id_u = PACK(box_id, touch_grid) ! hsml
            hsml_u = PACK(hsml, touch_grid) ! hsml
            mass_u = PACK(mass, touch_grid) ! masa
            dens_u = PACK(dens, touch_grid) ! densidad
            vol_u = mass_u / dens_u
            ! !vol_u = 4. * PI * hsml_u * hsml_u * hsml_u / 3. ! volumen

            do j = 1, 3
                pos_norm_u(:,j) = PACK(pos_norm(:,j), touch_grid) ! posición normalizada
                box_ijk_u(:,j) = PACK(box_ijk(:,j), touch_grid) ! (i,j,k) del box. De 0 a ng - 1
            end do
            
            ! Debido a que Aj aparece siempre multiplicado a vol, lo hacemos ahora
            !  y ya queda para el resto del cómputo.
            do j = 1, ninten
                Aj(j,:) = PACK(extra(:,j), touch_grid) * vol_u ! cantidad "extra" para intrínsecos
            end do
            do j = ninten + 1, n_extra
                Aj(j,:) = PACK(extra(:,j), touch_grid) ! cantidad "extra" para extrínsecos
            end do
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            ! PRINT*, "Contando la cantidad de partículas en cada box..."
            ! npart_per_cell = 0
            ! do i = 1, n_useful
            !     bin = box_id(i) + 1
            !     npart_per_cell(bin) = npart_per_cell(bin) + 1
            ! end do
            ! time_2 = OMP_GET_WTIME()
            ! PRINT '(2X,"Listo.")'
            ! PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            ! time_1 = OMP_GET_WTIME()

            ! PRINT*, "Generando array de cantidad de partículas acumualdas (previas) en boxes..."
            ! acum_part = 0
            ! do i = 2, ngrid3
            !     acum_part(i) = acum_part(i-1) + npart_per_cell(i-1)
            ! end do
            ! time_2 = OMP_GET_WTIME()
            ! PRINT '(2X,"Listo.")'
            ! PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            ! time_1 = OMP_GET_WTIME()

            ! PRINT*, "Generando array de índices (ordenados) de partículas en boxes..."
            ! id_part_in_cell = -1
            ! do i = 1, n_useful
            !     bin = box_id(i) + 1
            !     prev = acum_part(bin)
            !     do j = prev + 1, n_useful
            !         if (id_part_in_cell(j) .eq. -1) then
            !             exit
            !         end if
            !     end do
            !     id_part_in_cell(j) = i
            ! end do
            ! time_2 = OMP_GET_WTIME()
            ! PRINT '(2X,"Listo.")'
            ! PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            ! time_1 = OMP_GET_WTIME()

            PRINT*, "Calculando box amigo de cada partícula..."
            ! Calculamos a cuáles box toca cada partícula, según su hsml y posición
            do j = 1, 3
                do i = 1, n_useful
                    d2prev = MOD(pos_norm_u(i,j), 1.) ! que tan metida está en su box
                    if (d2prev < 0) then
                        d2prev = 1 + d2prev
                    end if
                    d2next = 1. - d2prev ! cuanto le falta para el box siguiente
                    first(i,j) = box_ijk_u(i,j) - FLOOR(ABS(1 + hsml_u(i) / h(j) - d2prev), kind=4) ! primero
                    last(i,j) = box_ijk_u(i,j) + FLOOR(ABS(1 + hsml_u(i) / h(j) - d2next), kind=4) ! último
                end do
            end do
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Calculando pesos Sj..."
            ! Calculamos los pesos Sj de cada partícula
            Sj = 0.
            !$OMP PARALLEL PRIVATE (i,bk,bj,bi,dx2,dy2,dz2,d,W)&
            !$OMP SHARED (first,last,pos_norm_u,hsml_u,Sj,hx,hy,hz)
            !$OMP DO SCHEDULE (STATIC)
                do i = 1, n_useful
                    do bk = first(i,3), last(i,3) ! del primero al último
                        dz2 = (pos_norm_u(i,3) - bk - 0.5) * hz
                        dz2 = dz2 * dz2
                        do bj = first(i,2), last(i,2) ! del primero al último
                            dy2 = (pos_norm_u(i,2) - bj - 0.5) * hy
                            dy2 = dy2 * dy2
                            do bi = first(i,1), last(i,1) ! del primero al último
                                dx2 = (pos_norm_u(i,1) - bi - 0.5) * hx
                                dx2 = dx2 * dx2
                                d = SQRT(dx2 + dy2 + dz2) ! desde su posición al centro del otro box
                                Sj(i) = Sj(i) + spline_d(d, hsml_u(i)) ! sumamos contribución
                            end do
                        end do
                    end do
                end do
            !$OMP END DO
            !$OMP END PARALLEL
            ! !Sj = Sj * vol_k !! Está comentado porque se hace maualmente a Ak al final (search: "Corrigiendo Pesos Ak")
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, " Escribiendo archivo 'Sj' con los pesos Sj de cada partícula...."
            OPEN(70, file="Sj", status="replace", form="unformatted", access="stream")
            WRITE(70) Sj
            CLOSE(70)
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Arreglando boxes vecinos..."
            ! Arreglamos First y Last
            do j = 1, 3
                do i = 1, n_useful
                    first(i,j) = MAX(MIN(first(i,j), ng(j) - 1), 0) ! primero
                    last(i,j) = MIN(MAX(last(i,j), 0), ng(j) - 1) ! último
                end do
            end do
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Calculando pesos Ak..."
            ! Calculamos los pesos Ak de cada box. Acá hay que tener en cuenta el valor de Aj ("extra") de cada partícula
            ! También calculamos las constantes de normalización Ik
            do i = 1, n_useful
                Aj_Sj(:,i) = Aj(:,i) / Sj(i) ! No molesta donde negativo, o Sj=0, porque no los usamos
            end do
            Ak = 0.
            Ik = 0.
            
            !$OMP PARALLEL PRIVATE (i,bk,bj,bi,bin,dx2,dy2,dz2,d,W)&
            !$OMP SHARED (box_id_u,Sj,first,last,pos_norm_u,hsml_u,Ak,Ik,Aj,Aj_Sj,vol_u,ngx,ngy,ngz,h2)
            !$OMP DO SCHEDULE (STATIC)
                do i = 1, n_useful
                    if (box_id_u(i) >= 0) then
                        if (Sj(i) > 0.) then ! en este caso debemos recorrer los boxes amigos
                            do bk = first(i,3), last(i,3) ! del primero al último
                                dz2 = (pos_norm_u(i,3) - bk - 0.5) * hz
                                dz2 = dz2 * dz2
                                do bj = first(i,2), last(i,2) ! del primero al último
                                    dy2 = (pos_norm_u(i,2) - bj - 0.5) * hy
                                    dy2 = dy2 * dy2
                                    do bi = first(i,1), last(i,1) ! del primero al último
                                        dx2 = (pos_norm_u(i,1) - bi - 0.5) * hx
                                        dx2 = dx2 * dx2
                                        d = SQRT(dx2 + dy2 + dz2) ! desde su posición al centro del otro box
                                        bin = ijk2boxid(bi, bj, bk, ngx, ngy, ngz) + 1
                                        W = spline_d(d, hsml_u(i))
                                        Ak(:,bin) = Ak(:,bin) + (W * Aj_Sj(:,i))
                                        Ik(bin) = Ik(bin) + (vol_u(i) * W) ! cte de normalización
                                    end do
                                end do
                            end do
                        else ! en este caso, solo consideramos el box más cercano (el suyo)
                            d = NORM2((pos_norm_u(i,:) - (/box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3)/) - 0.5) * h)
                            W = spline_d(d, hsml_u(i))
                            bin = ijk2boxid(box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3), ngx, ngy, ngz) + 1
                            Ak(:,bin) = Ak(:,bin) + Aj(:,i)
                            Ik(bin) = Ik(bin) + (vol_u(i) * W) ! cte de normalización
                        end if
                    else
                        if (Sj(i) > 0.) then ! en este caso debemos recorrer los boxes amigos
                            do bk = first(i,3), last(i,3) ! del primero al último
                                dz2 = (pos_norm_u(i,3) - bk - 0.5) * hz
                                dz2 = dz2 * dz2
                                do bj = first(i,2), last(i,2) ! del primero al último
                                    dy2 = (pos_norm_u(i,2) - bj - 0.5) * hy
                                    dy2 = dy2 * dy2
                                    do bi = first(i,1), last(i,1) ! del primero al último
                                        dx2 = (pos_norm_u(i,1) - bi - 0.5) * hx
                                        dx2 = dx2 * dx2
                                        d = SQRT(dx2 + dy2 + dz2) ! desde su posición al centro del otro box
                                        bin = ijk2boxid(bi, bj, bk, ngx, ngy, ngz) + 1
                                        W = spline_d(d, hsml_u(i))
                                        Ak(:,bin) = Ak(:,bin) + (W * Aj_Sj(:,i))
                                        Ik(bin) = Ik(bin) + (vol_u(i) * W) ! cte de normalización
                                    end do
                                end do
                            end do
                        end if
                    end if
                end do
            !$OMP END DO
            !$OMP END PARALLEL
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Corrigiendo pesos Ak..."
            do j = 1, ninten
                Ak(j,:) = Ak(j,:) / vol_k
            end do
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Normalizando valores de Ak según Ik..."
            ! Normalizamos según Ik
            if (ninten > 0) then
                do i = 1, ngrid3
                    if (Ik(i) > 0) then
                        do j = 1, ninten
                            Ak(j,i) = Ak(j,i) / Ik(i)
                        end do
                    end if
                end do
            end if
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Escribiendo archivo 'Ik' con pesos de normalización Ik..."
            OPEN(70, file="Ik", status="replace", form="unformatted", access="stream")
            WRITE(70) Ik
            CLOSE(70)
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            time_1 = OMP_GET_WTIME()
            
            PRINT*, "Dealocatando..."
            ! Dealocatamos
            deallocate(pos_norm_u, mass_u, vol_u, hsml_u)
            deallocate(box_ijk_u)
            deallocate(Sj, Aj, Aj_Sj)
            deallocate(first, last)
            ! Obtemos resultado final
            allocate(boxed_extra(n_extra,ngrid3))
            boxed_extra = Ak
            deallocate(Ak, Ik)
            time_2 = OMP_GET_WTIME()
            PRINT '(2X,"Listo.")'
            PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4," [s].")', time_2 - time_1 , time_2 - start_time
            PRINT*, ""
            PRINT*, "------ Fin subroutina 'grid_extra' ------"
            PRINT*, ""
        end subroutine grid_extra
    
end module box
