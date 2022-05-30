module box
    use kern
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
            end function boxid2k !0 -> ngrid³ -1
            
            integer*4 function boxid2j(id, nx, ny) result (j)
                !(box_id, nx, ny) -> j
                implicit none
                integer*4, intent(in) :: id, nx, ny

                j = MOD(id / nx, ny)
            end function boxid2j
            
            integer*4 function boxid2i(id, nx) result (i)
                !(box_id, nx) -> i
                implicit none
                integer*4, intent(in) :: id, nx

                i = MOD(id, nx)
            end function boxid2i
            
            subroutine boxid2ijk (id, nx, ny, i, j, k)
                !from (box_id, nx, ny) -> set: (i, j, k)
                implicit none
                integer*4, intent(in) :: id, nx, ny
                integer*4, intent(out) :: i, j, k

                i = boxid2i(id, nx)
                j = boxid2j(id, nx, ny)
                k = boxid2k(id, nx, ny)
            end subroutine boxid2ijk
            
            subroutine make_grid (nx, ny, nz, sx, sy, sz, center, grid)
                implicit none
                integer*4, intent(in) :: nx, ny, nz
                real*4, intent(in) :: sx, sy, sz
                real*4, dimension(3), intent(in) :: center
                real*4, dimension(:,:), allocatable, intent(out) :: grid
                integer*4 :: i, j, k, m
                real*4 :: rx, ry , rz
                
                rx = sx / (nx * 1.)
                ry = sy / (ny * 1.)
                rz = sz / (nz * 1.)
                allocate(grid(nx * ny * nz, 3))
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            m = ijk2boxid(i-1, j-1, k-1, nx, ny, nz) + 1
                            grid(m,1) = i * rx
                            grid(m,2) = j * ry
                            grid(m,3) = k * rz
                        end do
                    end do
                end do
                do m = 1, nx * ny * nz
                    grid(m,:) = grid(m,:) - 0.5 * (/sx + rx, sy + ry, sz + rz/) + center
                end do
            end subroutine make_grid
                    
            subroutine grid_extra (pos, mass, dens, hsml, extra, ninten, ngx, ngy, ngz, sx, sy, sz, center, boxed_extra)
                implicit none
                real*4, dimension(:), intent(in) :: mass
                real*4, dimension(SIZE(mass)), intent(in) :: dens, hsml
                real*4, dimension(SIZE(mass),3), intent(in) :: pos
                real*4, dimension(:, :), intent(in) :: extra ! (npart, n_extra)
                real*4, intent(in) :: sx, sy, sz
                integer*4, intent(in) :: ngx, ngy, ngz, ninten
                real*4, dimension(3), intent(in) :: center
                real*4, dimension(:,:), allocatable, intent(out) :: boxed_extra
                
                real*4, dimension(SIZE(mass),3) :: centred, pos_norm
                integer*4, dimension(SIZE(mass),3) :: box_ijk
                integer*4, dimension(SIZE(mass)) :: box_id
                logical, dimension(SIZE(mass)) :: in_grid, touch_grid

                real*4, dimension(3) :: h, s, s2, h2
                real*4 :: start_time, time_1, time_2
                integer*4, dimension(3) :: ng   
                real*4 :: hx, hy, hz, sx2, sy2, sz2, d2prev, d2next, d, W, vol_k
                integer*4 :: i, j, n_useful, ngrid3, n_extra, nexten
                integer*4 :: bi, bj, bk, bin

                real*4, allocatable, dimension(:) :: mass_u, dens_u, hsml_u, vol_u
                real*4, allocatable, dimension(:,:) :: pos_norm_u
                integer*4, allocatable, dimension(:,:) :: box_ijk_u, first, last
                real*4, allocatable, dimension(:) :: Sj, Ik
                real*4, allocatable, dimension(:,:) :: Aj ! este va a ser (n_extra, npart)
                real*4, allocatable, dimension(:,:) :: Ak ! este va a ser (n_extra, ng³)
                
                ! Checkeamos que "extra" esté bien
                if (SIZE(extra, 1) .ne. SIZE(mass)) then
                    PRINT*,     "Bad 'extra' input shape:"
                    PRINT*,     " ", SIZE(extra, 1), " was given, and ", SIZE(mass), " was expected"
                    PRINT*,     "Exiting."
                    call EXIT(1)
                end if
                if (SIZE(extra, 2) .lt. ninten) then
                    PRINT*,     "Bad 'extra' input shape:"
                    PRINT*,     " ", SIZE(extra, 2), " was given, and not lower than", ninten, " was expected"
                    PRINT*,     "Exiting."
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
                
                call CPU_TIME (start_time)
                call CPU_TIME (time_1)
                
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
                
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                PRINT*, ""
                call CPU_TIME (time_1)
                
                PRINT*, "Calculando índice de box de partículas..."
                ! Pasamos de (i,j,k) a box_id
                do i = 1, SIZE(mass) 
                    if (in_grid(i)) then ! las que están dentro de la grilla...
                        box_id(i) = ijk2boxid(box_ijk(i, 1), box_ijk(i, 2), box_ijk(i, 3), ngx, ngy, ngz) !va de 0 a ng-1
                    else 
                        box_id(i) = -1 !-1 tienen las que caen fuera de la grilla
                    end if
                end do
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, "Escribiendo archivo 'box_id' con índice de box de partículas..."
                OPEN(90, file="box_id", status="replace", form="unformatted", access="stream")
                WRITE(90) box_id
                CLOSE(90)
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)

                PRINT*, "Reduciendo set de datos..."
                ! Cantidad de partículas que caen dentro de la grilla
!                 touch_grid = in_grid
                n_useful = COUNT(touch_grid) ! esta es la cantidad de partículas ("útiles") que usaremos
                PRINT '(3X,"Partículas dentro del grid           : ",I10)', COUNT(in_grid)
                PRINT '(3X,"Dentro/Total                         : ",F6.4)', COUNT(in_grid) / (1. * SIZE(mass))
                PRINT '(3X,"Partículas interactuantes con el grid: ",I10)', n_useful
                PRINT '(3X,"Dentro/Interactuantes                : ",F6.4)', COUNT(in_grid) / (1. * n_useful)
                PRINT '(3X,"Interactuantes/Total                 : ",F6.4)', n_useful / (1. * SIZE(mass))
                
                ! Alocatamos 
                allocate(mass_u(n_useful), dens_u(n_useful), hsml_u(n_useful), vol_u(n_useful)) 
                allocate(pos_norm_u(n_useful,3), box_ijk_u(n_useful,3))
                allocate(first(n_useful,3), last(n_useful,3))
                allocate(Sj(n_useful))
                allocate(Aj(n_extra,n_useful))
                allocate(Ak(n_extra,ngrid3), Ik(ngrid3))
                
                ! Ponemos en estos array los datos correspondientes.
                hsml_u = PACK(hsml, touch_grid) ! hsml
                mass_u = PACK(mass, touch_grid) ! masa
                dens_u = PACK(dens, touch_grid) ! densidad
                vol_u = mass_u / dens_u
!                 vol_u = 4. * PI * hsml_u * hsml_u * hsml_u / 3. ! volumen

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
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
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
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, "Calculando pesos Sj..."
                ! Calculamos los pesos Sj de cada partícula
                Sj = 0.
                do i = 1, n_useful
                    do bk = first(i,3), last(i,3) ! del primero al último
                        do bj = first(i,2), last(i,2) ! del primero al último
                            do bi = first(i,1), last(i,1) ! del primero al último
                                d = SQRT(SUM((pos_norm_u(i,:) - (/bi, bj, bk/) - 0.5)**2 * h2)) ! desde su posición al centro del otro box
                                Sj(i) = Sj(i) + spline_d(d, hsml_u(i)) ! sumamos contribución
                            end do
                        end do
                    end do
                end do
!                 Sj = Sj * vol_k !! Está comentado porque se hace maualmente a Ak al final
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, " Escribiendo archivo 'Sj' con los pesos Sj de cada partícula...."
                OPEN(70, file="Sj", status="replace", form="unformatted", access="stream")
                WRITE(70) Sj
                CLOSE(70)
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, "Arreglando boxes vecinos..."
                ! Arreglamos First y Last
                do j = 1, 3
                    do i = 1, n_useful
                        first(i,j) = MAX(MIN(first(i,j), ng(j) - 1), 0) ! primero
                        last(i,j) = MIN(MAX(last(i,j), 0), ng(j) - 1) ! último
                    end do
                end do
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, "Calculando pesos Ak..."
                ! Calculamos los pesos Ak de cada box. Acá hay que tener en cuenta el valor de Aj ("extra") de cada partícula
                ! También calculamos las constantes de normalización Ik
                do i = 1, n_useful
                    if (box_id(i) >= 0) then
                        if (Sj(i) > 0.) then ! en este caso debemos recorrer los boxes amigos
                            do bk = first(i,3), last(i,3) ! del primero al último
                                do bj = first(i,2), last(i,2) ! del primero al último
                                    do bi = first(i,1), last(i,1) ! del primero al último
                                        bin = ijk2boxid(bi, bj, bk, ngx, ngy, ngz) + 1
                                        d = SQRT(SUM((pos_norm_u(i,:) - (/bi, bj, bk/) - 0.5)**2 * h2))
                                        W = spline_d(d, hsml_u(i))
                                        Ak(:,bin) = Ak(:,bin) + (W * Aj(:,i) / Sj(i))
                                        Ik(bin) = Ik(bin) + (vol_u(i) * W) ! cte de normalización
                                    end do
                                end do
                            end do
                        else ! en este caso, solo consideramos el box más cercano (el suyo)
                            d = SQRT(SUM((pos_norm_u(i,:) - (/box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3)/) - 0.5)**2 * h2))
                            W = spline_d(d, hsml_u(i))
                            bin = ijk2boxid(box_ijk_u(i,1), box_ijk_u(i,2), box_ijk_u(i,3), ngx, ngy, ngz) + 1
                            Ak(:,bin) = Ak(:,bin) + Aj(:,i)
                            Ik(bin) = Ik(bin) + (vol_u(i) * W) ! cte de normalización
                        end if
                    else
                       if (Sj(i) > 0.) then ! en este caso debemos recorrer los boxes amigos
                            do bk = first(i,3), last(i,3) ! del primero al último
                                do bj = first(i,2), last(i,2) ! del primero al último
                                    do bi = first(i,1), last(i,1) ! del primero al último
                                        bin = ijk2boxid(bi, bj, bk, ngx, ngy, ngz) + 1
                                        d = SQRT(SUM((pos_norm_u(i,:) - (/bi, bj, bk/) - 0.5)**2 * h2))
                                        W = spline_d(d, hsml_u(i))
                                        Ak(:,bin) = Ak(:,bin) + (W * Aj(:,i) / Sj(i))
                                        Ik(bin) = Ik(bin) + (vol_u(i) * W) ! cte de normalización
                                    end do
                                end do
                            end do
                        end if  
                    end if
                end do
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, "Corrigiendo pesos Ak..."
                do j = 1, ninten
                    Ak(j,:) = Ak(j,:) / vol_k
                end do
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
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
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, "Escribiendo archivo 'Ik' con pesos de normalización Ik..."
                OPEN(70, file="Ik", status="replace", form="unformatted", access="stream")
                WRITE(70) Ik
                CLOSE(70)
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                call CPU_TIME (time_1)
                
                PRINT*, "Dealocatando..."
                ! Dealocatamos
                deallocate(pos_norm_u, mass_u, vol_u, hsml_u)
                deallocate(box_ijk_u)
                deallocate(Sj, Aj)
                deallocate(first, last)
                ! Obtemos resultado final
                allocate(boxed_extra(n_extra,ngrid3))
                boxed_extra = Ak
                deallocate(Ak, Ik)
                call CPU_TIME (time_2)
                PRINT '(2X,"Listo.")'
                PRINT '(2X,"Este bloque: ",F8.4," [s]. Total: ",F8.4)', time_2 - time_1 , time_2 - start_time
                PRINT*, ""
                PRINT*, "------ Fin subroutina 'grid_extra' ------"
                PRINT*, ""
            end subroutine grid_extra
    
end module box
