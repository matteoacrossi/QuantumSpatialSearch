function apply_dyson!(H, psi, dt, dysonOrder)
    kq = copy(psi)
    # evolution (dysonOrder is the maximum order of dt)
    for k = 1 : dysonOrder
        kq = - 1im * dt/k * H * kq
        psi .+= kq
    end
end
