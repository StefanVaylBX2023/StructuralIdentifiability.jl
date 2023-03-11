

function integral_tester(ode::ODE{P}) where P <: MPolyElem
    ord = 10
    ys = ode.y_vars
    result = perform_substitution(ode)
    init_ode = result[3]
    res_ode = result[1]
    integrals = result[2]
    map_x = Dict(x => fmpq(rand(1:10)) for x in init_ode.x_vars)
    map_p = Dict(p => fmpq(rand(1:10)) for p in init_ode.parameters)
    map_u = Dict(u => [fmpq(rand(1:100)) for i in 0:ord] for u in init_ode.u_vars)
    ps_ode = power_series_solution(init_ode, map_p, map_x, map_u, ord)
    #----------------------------------------------------------------------------------------------------------------------------------------
    map_x_str = Dict(var_to_str(x) => f for (x,f) in map_x)
    map_u_str = Dict(var_to_str(u) => f for (u,f) in map_u)
    map_p_str = Dict(var_to_str(p) => f for (p,f) in map_p)
    #----------------------------------------------------------------------------------------------------------------------------------------
    cvars = Dict{fmpq_mpoly, fmpq}()
    for (x,f) in integrals
        var = vars(f)
        temp= Dict{fmpq_mpoly, fmpq_mpoly}()
        for v in var
            t = var_to_str(v)
            temp[v] = map_x_str[t]*one(res_ode.poly_ring)
        end
        cvars[x] = first(coefficients(eval_at_dict(f, temp)))
    end
    map_p_res = Dict{fmpq_mpoly, fmpq}(p => (map_p_str[var_to_str(p)]) for p in res_ode.parameters if var_to_str(p) ∉ map(var_to_str, collect(keys(cvars))))
    for (x, f) in cvars
        map_p_res[parent_ring_change(x, res_ode.poly_ring)] = f
    end
    map_u_res = Dict(u => (map_u_str[var_to_str(u)]) for u in res_ode.u_vars)
    map_x_res = Dict(x => (map_x_str[var_to_str(x)]) for x in res_ode.x_vars) 
    res_ps_ode = power_series_solution(res_ode, map_p_res, map_x_res, map_u_res, ord)
    #----------------------------------------------------------------------------------------------------------------------------------------
    result = false
    for y in ys
        y = var_to_str(y)
        init_y = str_to_var(y, init_ode.poly_ring)
        res_y = str_to_var(y, res_ode.poly_ring)
        result = (ps_ode[init_y] == res_ps_ode[res_y])
        if result == false 
            return false
        end
    end
    return true
end

@testset "test integrals" begin
    cases = [@ODEmodel(
        x1'(t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k4 * x6(t),
        x2'(t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k3 * x4(t),
        x3'(t) = k3 * x4(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
        x4'(t) = k1 * x1(t) * x2(t) - k2 * x4(t) - k3 * x4(t),
        x5'(t) = k4 * x6(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
        x6'(t) = -k4 * x6(t) - k5 * x6(t) + k6 * x3(t) * x5(t),
        y1(t) = x3(t)
    ), @ODEmodel(
        x1'(t) = -x1(t)//x2(t) + x2(t)//x1(t),
        x2'(t) = -x2(t)//x1(t) + x1(t)//x2(t),
        y1(t) = x1(t),
        y2(t) = x2(t)
    ),@ODEmodel(
        x1'(t) = (-x1(t)-x2(t))//x4(t),
        x2'(t) = x2(t)//x4(t),
        x3'(t) = x1(t)//x4(t),
        x4'(t) = x4(t),
        y1(t) = x1(t)
    ),@ODEmodel(
        S'(t) = -b * ln(t) * S(t) // (S(t) + L(t) + ln(t) + R(t) + Q(t)) - u(t) * S(t) // (S(t) + L(t) + ln(t) + R(t) + Q(t)),
        L'(t) = b * ln(t) * S(t) // (S(t) + L(t) + ln(t) + R(t) + Q(t)) - a * L(t),
        ln'(t) = a * L(t) - g * ln(t) + s * Q(t),
        R'(t) = u(t) * S(t) // (S(t) + L(t) + ln(t) + R(t) + Q(t)) + e * g * ln(t),
        Q'(t) = (1 - e) * g * ln(t) - s * Q(t),
        y1(t) = ln(t) // (S(t) + L(t) + ln(t) + R(t) + Q(t))
    ),@ODEmodel(
        Sh'(t) = mh * (Sh(t) + Lh(t) + Ih(t) + Rh(t)) - bv * Sh(t) * Iv(t) // (Sh(t) + Lh(t) + Ih(t) + Rh(t)) - mh * Sh(t),
        Lh'(t) = bv * Sh(t) * Iv(t) // (Sh(t) + Lh(t) + Ih(t) + Rh(t)) - ah * Lh(t) - mh * Lh(t) ,
        Ih'(t) = ah * Lh(t) - gh * Ih(t) - mh * Ih(t),
        Rh'(t) = gh * Ih(t) - mh * Rh(t),
        Sv'(t) = mv * (Sv(t) + Lv(t) + Iv(t))- bh * Sv(t) * Ih(t) // (Sh(t) + Lh(t) + Ih(t) + Rh(t)) - mv * Sv(t),
        Lv'(t) = bh * Sv(t) * Ih(t) // (Sh(t) + Lh(t) + Ih(t) + Rh(t)) - av * Lv(t) - mv * Lv(t),
        Iv'(t) = av * Lv(t) - mv * Iv(t),
        y1(t) = Ih(t) // (Sh(t) + Lh(t) + Ih(t) + Rh(t))
    
    ),@ODEmodel(
        x1'(t) = -x1(t) + x2(t),
        x2'(t) = x1(t),
        y1(t) = x1(t) - x2(t)
    ),@ODEmodel(
    KS00'(t) = -a00 * K(t) * S00(t) + b00 * KS00(t) + gamma0100 * FS01(t) + gamma1000 * FS10(t) + gamma1100 * FS11(t),
    KS01'(t) = -a01 * K(t) * S01(t) + b01 * KS01(t) + c0001 * KS00(t) - alpha01 * F(t) * S01(t) + beta01 * FS01(t) + gamma1101 * FS11(t),
    KS10'(t) = -a10 * K(t) * S10(t) + b10 * KS10(t) + c0010 * KS00(t) - alpha10 * F(t) * S10(t) + beta10 * FS10(t) + gamma1110 * FS11(t),
    FS01'(t) = -alpha11 * F(t) * S11(t) + beta11 * FS11(t) + c0111 * KS01(t) + c1011 * KS10(t) + c0011 * KS00(t),
    FS10'(t) = a00 * K(t) * S00(t) - (b00 + c0001 + c0010 + c0011) * KS00(t),
    FS11'(t) = a01 * K(t) * S01(t) - (b01 + c0111) * KS01(t),
    K'(t) = a10 * K(t) * S10(t) - (b10 + c1011) * KS10(t),
    F'(t) = alpha01 * F(t) * S01(t) - (beta01 + gamma0100) * FS01(t),
    S00'(t) = alpha10 * F(t) * S10(t) - (beta10 + gamma1000) * FS10(t),
    S01'(t) = alpha11 * F(t) * S11(t) - (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
    S10'(t) = -a00 * K(t) * S00(t) + (b00 + c0001 + c0010 + c0011) * KS00(t) - a01 * K(t) * S01(t) + (b01 + c0111) * KS01(t) - a10 * K(t) * S10(t) + (b10 + c1011) * KS10(t),
    S11'(t) = -alpha01 * F(t) * S01(t) + (beta01 + gamma0100) * FS01(t) - alpha10 * F(t) * S10(t) + (beta10 + gamma1000) * FS10(t) - alpha11 * F(t) * S11(t) + (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
    y1(t) = F(t),
    y2(t) = S00(t),
    y3(t) = S01(t),
    y4(t) = S10(t),
    y5(t) = S11(t)
    )
    ]

    for c in cases
        
        @test (integral_tester(c) == true)
    end
end
