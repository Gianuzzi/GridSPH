module kern

    implicit none
    real*8, parameter :: PI = 4. * ATAN(1.)
    real*8, parameter :: ONE_PI = 1. / PI

     contains
     
        real*8 function get_q(x, y, z, h) result (q)
            implicit none
            real*8, intent(in) :: x, y, z, h
            q = SQRT(x*x + y*y + z*z) / h
        end function get_q
        
        real*8 function spline_q(q, h) result (w)
            implicit none
            real*8, intent(in) :: q, h
            real*8 :: aux
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
        
        real*8 function spline_d(d, h) result (w)
            implicit none
            real*8, intent(in) :: d, h
            w = spline_q(d / h, h)
        end function spline_d
        
        real*8 function spline_xyz(x, y, z, h) result (w)
            implicit none
            real*8, intent(in) :: x, y, z, h
            w = spline_q(get_q(x, y, z, h), h)
        end function spline_xyz
    
end module kern
