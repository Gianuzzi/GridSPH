module kern

    implicit none
    real, parameter :: PI = 4. * ATAN(1.)
    real, parameter :: ONE_PI = 1. / PI

     contains
     
        real function get_q(x, y, z, h) result (q)
            implicit none
            real, intent(in) :: x, y, z, h
            q = SQRT(x*x + y*y + z*z) / h
        end function get_q
        
        real function spline_q(q, h) result (w)
            implicit none
            real, intent(in) :: q, h
            real :: aux
            if (q <= 1) then
                w = 1 - 1.5 * q * q * (1 - q * 0.5)
            else if ((1 < q) .and. (q <= 2)) then
                aux = (2 - q)
                w = 0.25 * aux * aux * aux
            else
                w = 0
            end if
            w = w * ONE_PI / (h * h * h)
        end function spline_q
        
        real function spline_d(d, h) result (w)
            implicit none
            real, intent(in) :: d, h
            w = spline_q(d / h, h)
        end function spline_d
        
        real function spline_xyz(x, y, z, h) result (w)
            implicit none
            real, intent(in) :: x, y, z, h
            w = spline_q(get_q(x, y, z, h), h)
        end function spline_xyz
    
end module kern
