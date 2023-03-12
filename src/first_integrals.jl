


# -----------------------------------------------------------------------------


"""
    coeff_matrix(xeq)
Takes x_equations as input and returns matrix of coefficients of allmonomials in x_equations.

Input:
- `xeq` - x_equations

Output: 
-  fmpq_mat object containing coefficients of all monomials in x_equations and mappings from equations to rows in matrix
"""


function coeff_matrix(xeq::Dict{fmpq_mpoly, <:Union{Generic.Frac{fmpq_mpoly}, fmpq_mpoly}})
    eq2denom = Dict{fmpq_mpoly, fmpq_mpoly}(x => unpack_fraction(f)[2] for (x, f) in xeq)
    denom = lcm(collect(values(eq2denom)))
    ode_x = Dict{fmpq_mpoly, fmpq_mpoly}([x => unpack_fraction(f * denom)[1] for (x, f) in xeq])
    equations = collect(values(ode_x))
    monoms = Set{fmpq_mpoly}(reduce(union, map(monomials, equations)))
    matrix = zero(Nemo.MatrixSpace(Nemo.QQ, length(xeq), length(monoms)))
    map_monom = Dict{fmpq_mpoly, Int64}(monomial => i for (i, monomial) in enumerate(monoms))
    map_eq = Dict{fmpq_mpoly, Int64}(x => i for (i, (x, f)) in enumerate(ode_x))
    for (x, f) in ode_x
        for (c, m) in zip(coefficients(f), monomials(f))
            matrix[map_eq[x], map_monom[m]] = c
        end
    end
    return matrix, map_eq
end

# -----------------------------------------------------------------------------

"""
    find_first_integrals(ode, out)
This function takes ode system as an input and by default returns
an array of first integrals solutions available at the moment.
If one specifies out = "map" it returns a dictionary where keys are variables and values are coefficients in the first integral solution.

Input:
- `ode` - ode system
- 'out' - default value is "polynom" 
          optional  value is "map" 
Output: 
-  Array of first integral solutions or dictionary with keys as variables and values as coefficients in the first integral solution.
"""

function find_first_integrals(ode::ODE{P}, out = :polynom) where P <: MPolyElem
    xeq = ode.x_equations
    matrix, map_eq = coeff_matrix(xeq)
    n_col, solution_set = Nemo.nullspace(transpose(matrix))
    (n_col == 0) && return []
    solution_matrix = transpose(solution_set)[1, : ]
    solution_matrices = [transpose(solution_set)[i, : ] for i in 1:n_col]
    solution = Dict{fmpq_mpoly,fmpq}(Dict(x => solution_matrix[i] for (x, i) in map_eq if solution_matrix[i] != 0))
    variables = keys(solution)
    ring = parent(first(variables))
    sol = ring(first(variables) * solution[first(variables)])
    sol_poly = Array{fmpq_mpoly, 1}()
    for sol in solution_matrices
        dict_sol = Dict{fmpq_mpoly, fmpq}(Dict(x => sol[i] for (x, i) in map_eq if sol[i] != 0))
        poly = ring(first(keys(dict_sol)) * dict_sol[first(keys(dict_sol))])
        for (x, coef) in dict_sol
            if x != first(variables)
                poly += ring(x * coef)
            end
        end
        push!(sol_poly, poly)
    end
    if out == :map
        return solution
    elseif out == :polynom
        return sol_poly
    else 
        throw("Wrong output format $out.")
    end
end




# -----------------------------------------------------------------------------

"""
    get_max_exp(xeq,var)
Finds maximum exponent of variable `var` in all x_equations.

"""

function get_max_exp(xeq::Dict{fmpq_mpoly, <: Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}, var::fmpq_mpoly)
    xs = keys(xeq)
    map_exp = Dict{fmpq_mpoly, Int64}(x => 0 for x in xs)
    for (x, f) in xeq
        f = unpack_fraction(f)
        for variable in xs
            map_exp[variable] = max(map_exp[variable], degree(f[1], variable), degree(f[2], variable))
        end
    end
    return map_exp[var]
end

# -----------------------------------------------------------------------------

"""
    get_occurences(xeq,var)
Finds occurences of variable in monomials, tottal occurences  and maximum coefficient of variable in the whole ODE system.

"""

function get_occurences(xeq::Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}, var::fmpq_mpoly)
    total_occur = 0
    monom_occur = 0
    max_coef = 0
    for (x, f) in xeq
        monoms = union(monomials(unpack_fraction(f)[1]), monomials(unpack_fraction(f)[2]))
        terms_ = union(terms(unpack_fraction(f)[1]), terms(unpack_fraction(f)[2]))
        for x in monoms
            if var in vars(x)
                total_occur += 1  
            end
            if (var in vars(x)) && (length(vars(x)) >= 2)
                monom_occur += 1
            end
        end
        for x in terms_
            if var in vars(x) && abs(collect(coefficients(x))[1]) > max_coef
                max_coef = abs(collect(coefficients(x))[1])
            end
        end
    end
    return monom_occur, total_occur, max_coef
end

# -----------------------------------------------------------------------------

function map_var2par(xeq::Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}, candidates::Vector{fmpq_mpoly}, scoring_functions)
    ranks = Dict(x => collect(Base.Iterators.flatten([F(xeq, x) for F in scoring_functions])) for x in candidates)
    return ranks
end

# -----------------------------------------------------------------------------

"""
    find_var2sub(mappings)
Based on results of solution finds variable which suits best for substitution.
"""

function find_var2sub(mappings)
    to_compare = collect(values(mappings))
    ([to_compare[i][1] for i in 1:lastindex(to_compare)] == [0 for i in 1:lastindex(to_compare)]) && return first(keys(mappings))

    to_compare = [to_compare[i] for i in 1:length(to_compare) if to_compare[i][1] > 0]
    (to_compare == []) && return 0
    
    (length(to_compare) == 1) && return findfirst(x -> x == to_compare[1], mappings)

    minimum = min(to_compare...)

    return findfirst(x -> x == minimum, mappings)
end

# -----------------------------------------------------------------------------

function get_var2sub(xeq::Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}, candidates::Vector{fmpq_mpoly}, scoring_functions)
    mappings = map_var2par(xeq, candidates, scoring_functions)
    return find_var2sub(mappings)
end

# -----------------------------------------------------------------------------

"""
    construct_substitution(ode)
Main function which proceeds with substitution based on first possible solution by constructing and applying 
substitution based on the solution.

Input:
- `ode` - initial ODE system
- `ode_aux` - auxiliary ODE system which construct another system 
based on actions applied to main one but without substitution of c_i parameters

Output: 
- `new_ode` - ODE system with applied substitution
- `integral` - Dictionary with current first integral where key is a paramter and value is a first integral
- `ode_aux` - auxiliary ODE system without substitution of c_i parameters
"""



function construct_substitution(ode::ODE{P}, ode_aux::ODE{P}) where P <: MPolyElem
    
    new_xeqs = y_integrals(ode)
    ode, new_xs = add_x_eqs(ode, new_xeqs)
    ode = add_y_integrals(ode, new_xs)

    new_xeqs_aux = y_integrals(ode_aux)
    ode_aux, new_xs_aux = add_x_eqs(ode_aux, new_xeqs_aux)
    ode_aux = add_y_integrals(ode_aux, new_xs_aux)
    
    new_const = gen_new_name(ode, "c")
    ode = add_parameter(ode, new_const)

    solution = find_first_integrals(ode, :map)  
    (solution == []) && return 0, [], 0

    candidates = collect(keys(solution))
    var2sub = get_var2sub(ode.x_equations, candidates, [get_max_exp, get_occurences])
    ring = ode.poly_ring

    new_const = ring(str_to_var(new_const, ring))
    integral = Dict()
    integral[new_const] = ring(solution[var2sub] * var2sub)
    lc = length(candidates)
    for i in 1:lc
        if candidates[i] != var2sub
            integral[new_const] = integral[new_const] + ring(candidates[i] * solution[candidates[i]])
        end
    end

    sub = new_const
    for i in 1:lc
        if candidates[i] != var2sub
            sub = sub - ring(candidates[i] * solution[candidates[i]])
        end
    end

    new_x = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}()
    new_y = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}()
    new_u = copy(ode.u_vars)

    sub = (1 // solution[var2sub]) * sub
    for (x, f) in ode.x_equations
        if x != var2sub
            new_x[x] = make_substitution(f, var2sub, sub)
        end
    end
    for (y,f) in ode.y_equations
        new_y[y] = make_substitution(f, var2sub, sub)
    end
    
    new_vars = map(var_to_str, setdiff(gens(ring), [var2sub]))
    S, _ = Nemo.PolynomialRing(Nemo.QQ, new_vars)
    dict_type = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}
    fin_x = dict_type(parent_ring_change(x, S) => parent_ring_change(f, S) for (x, f) in new_x)
    fin_y = dict_type(parent_ring_change(y, S) => parent_ring_change(f, S) for (y, f) in new_y)
    fin_u = Array{fmpq_mpoly,1}([parent_ring_change(u, S) for u in new_u])

    new_ode =  ODE{fmpq_mpoly}(fin_x, fin_y, fin_u)
    return new_ode, integral, ode_aux
end

# -----------------------------------------------------------------------------



"""
    perform_substitution(ode)
Main function which recurcively itterates through all the possible substitutions by applying function construct_substitution.

Input:
- `ode` - initial ODE system

Output: 
- final ODE system after all the substitutions have been applied
- `integrals` - Dictionary with all the first integrals where keys are paramters and values are first integrals
- `ode_aux` - auxiliary ODE system without substitution of c_i parameters
"""
function perform_substitution(ode::ODE{P}) where P <: fmpq_mpoly
    ys = collect(keys(ode.y_equations))
    integrals = Dict{P, Union{P, Generic.Frac{P}}}()
    ode = add_u_eqs(ode)
    ode_aux = ode
    subbed, integral, ode_aux = construct_substitution(ode,ode_aux)
    temp = subbed
    temp2 = ode_aux
    if subbed == 0
        @debug "There are no substitutions for this case"
        return ode, integrals, ode_aux
    else
        integrals[first(keys(integral))] = integral[first(keys(integral))]
        while subbed != 0
            temp = subbed
            temp2 = ode_aux
            subbed, integral, ode_aux = construct_substitution(subbed, ode_aux)
            if integral != []
                integrals[first(keys(integral))] = integral[first(keys(integral))]
            end
        end
    end
    return dfs_cleanup(temp, ys), integrals, temp2
end

# -----------------------------------------------------------------------------

"""
    add_x_eqs(ode, extra_x)
Adds for all y_equations new x_equations to involve them in substitution process.
"""

function add_x_eqs(ode::ODE{P}, extra_x) where P <: MPolyElem

    new_var_names =
        vcat(collect(map(var_to_str, gens(ode.poly_ring))), collect(keys(extra_x)))
    new_ring, _ = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_var_names)
    new_x_eqs = Dict{P, Union{P, Generic.Frac{P}}}(parent_ring_change(x, new_ring) => parent_ring_change(f, new_ring) for (x, f) in ode.x_equations)
    new_y_eqs = Dict{P, Union{P, Generic.Frac{P}}}(parent_ring_change(y, new_ring) => parent_ring_change(g, new_ring) for (y, g) in ode.y_equations)
    extra_x_eqs = Dict{P, Union{P, Generic.Frac{P}}}(str_to_var(x, new_ring) => parent_ring_change(g, new_ring) for (x, g) in extra_x)
    new_xs = collect(keys(extra_x_eqs))
    merge!(new_x_eqs, extra_x_eqs)
    new_us = map(v -> switch_ring(v, new_ring), ode.u_vars)
    return ODE{P}(new_x_eqs, new_y_eqs, new_us), new_xs
end

# -----------------------------------------------------------------------------

"""
    y_integrals(ode)
    Finds x_equations which are equal to some y's and could be used in substitution.
"""

function y_integrals(ode::ODE{P}) where P <: MPolyElem
    new_xeq = Dict{String,Union{P, Generic.Frac{P}}}()
    yeq = copy(ode.y_equations)
    xeq = copy(ode.x_equations)
    for (y,f) in yeq
        if f ∉ collect(values(xeq))
            new_x = gen_new_name(ode, "x")
            ode = add_parameter(ode, new_x)
            new_xeq[new_x] = f
        end
    end
    return new_xeq
end

# -----------------------------------------------------------------------------



"""
    add_y_integrals(ode, new_x)
For each new x, adds a new y and adds the equation y = x.

"""

function add_y_integrals(ode::ODE{P}, new_x::Array{fmpq_mpoly, 1}) where P <: MPolyElem
    ring = ode.poly_ring
    for x in new_x
        x = parent_ring_change(x, ring)
        new_y = gen_new_name(ode, "y")
        ode = add_outputs(ode, Dict(new_y => x))
    end
    return ode
end

# -----------------------------------------------------------------------------

"""
    add_u_eqs(ode)
For each u in ode.u_vars, adds a new x and y, and adds the equation x = 1.
"""


function add_u_eqs(ode::ODE{P}) where P <: MPolyElem

    us = copy(ode.u_vars)
    
    new_y = gen_new_name(ode, "y")
    new_x = gen_new_name(ode, "x")
    ode, _ = add_x_eqs(ode, Dict(new_x => one(ode.poly_ring)))
    ode = add_outputs(ode, Dict(new_y => str_to_var(new_x, ode.poly_ring)))
    
    (us == []) && return ode

    for u in us
        new_y = gen_new_name(ode, "y")
        new_x = gen_new_name(ode, "x")
        ode, _ = add_x_eqs(ode, Dict(new_x => u))
        ode = add_outputs(ode, Dict(new_y => str_to_var(new_x, ode.poly_ring)))
    end    

    return ode
end

# -----------------------------------------------------------------------------

"""
    gen_new_name(ode::ODE, var)
Generates new name for variable/constant not involved in the ODE
"""

function gen_new_name(ode::ODE{P}, var::String) where P <: MPolyElem
    str = map(var_to_str, gens(ode.poly_ring))
    counter = 1
    while ("$var$counter") ∈ str
        counter += 1
    end
    return "$var$counter" 
end

# -----------------------------------------------------------------------------

"""
    dfs_cleanup(ode, ys)

Cleans up the ODE system after the substitution process.

"""

function dfs_cleanup(ode::ODE{P}, ys::Vector{fmpq_mpoly}) where P <: MPolyElem
    ring = ode.poly_ring
    ys = [parent_ring_change(y, ring) for y in ys]
    graph = construct_graph(ode)
    raw = traverse_outputs(graph, ys)
    valid_vars = (union(collect(values(raw))...))
    vars_str = [var_to_str(x) for x in valid_vars]
    for (y, f) in ode.y_equations
        if y ∉ valid_vars
            temp = map(var_to_str, vars(f))
            if issubset(temp, vars_str)
                push!(valid_vars, y)
                push!(vars_str, var_to_str(y))
            end
        end
    end
    S, _ = Nemo.PolynomialRing(Nemo.QQ, vars_str)
    dict_type = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}
    new_xeq = dict_type(parent_ring_change(x, S) => parent_ring_change(f, S) for (x, f) in ode.x_equations if x in valid_vars)
    new_yeq = dict_type(parent_ring_change(y, S) => parent_ring_change(f, S) for (y, f) in ode.y_equations if y in valid_vars)
    new_u = Array{fmpq_mpoly, 1}([parent_ring_change(u, S) for u in ode.u_vars if u in valid_vars])
    new_ode = ODE{fmpq_mpoly}(new_xeq, new_yeq, new_u)
    return new_ode
end

