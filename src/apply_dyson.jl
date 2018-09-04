function apply_dyson!(H, psi, dt, dysonOrder)
    kq = copy(psi)
    tmp = similar(kq)
    # evolution (dysonOrder is the maximum order of dt)
    for k = 1 : dysonOrder
        #kq = - 1im * dt/k * H * kq
        mul!(tmp, H, kq, -1im * dt/k, 0)
        psi .+= tmp
        copyto!(kq, tmp)
    end
end
