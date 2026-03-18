


function flang_l2_norm(n,u,v) result(res)
    implicit none
    integer :: n
    real(kind(1.0d0)) :: u(n), v(n), res
    res = norm2(hypot(u, v))
end function