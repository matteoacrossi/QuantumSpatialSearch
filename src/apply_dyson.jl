function apply_dyson!(H, psi, dt::Real, dysonOrder::Integer)
    kq = copy(psi)
    tmp = similar(kq)
    for k = 1 : dysonOrder
        mul!(tmp, H, kq, -1im * dt/k, 0)
        psi .+= tmp
        copyto!(kq, tmp)
    end
end
