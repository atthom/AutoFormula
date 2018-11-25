using Base.MathConstants

import Base.==

# https://fr.wikipedia.org/wiki/Table_de_constantes_math%C3%A9matiques
cmm = 0.2614972128
cbern = 0.280169499
cGKW = 0.3036630029
cHSM = 0.35323637
comega = 0.5671432904
cEM = 0.5772156649
cCD = 0.6243299885
cC = 0.6434105462
cpj = 0.6601618158
crobb = 0.6617071822
clp = 0.6627434193
cET = 0.70258
cLR = 0.7642236535
cGauss = 0.834626842
cBrun = 0.8705883800
cK = 0.9159655941
cKL = 1.1865691104
cApery = 1.2020569031
cGK = 1.28242712
cConway = 1.30357729
cMills = 1.3063778838
cPlastique = 1.3247179572
csqrt2 = 1.41421356237
cRS = 1.4513692348
cEB = 1.6066951524
cGolden = 1.6180339887
cKhintchine = 2.685452001
cFeigenbaum = 4.6692016091 

math = [cRS, cEB, cGolden, csqrt2, cPlastique, cMills, cConway, cGK, 
    cApery, cKL, cK, cBrun, cGauss, cLR, cET, clp, crobb,
    cpj, cC, cCD, cEM, comega, cHSM, cGKW, cbern, cmm,
    cKhintchine, cFeigenbaum]


coulomb = 8.987551787e9
ce = 1.602176620e-19
cG = 6.67408e-11
cpermeability_magne = 1.2566370614e-6
cpermeability_elec = 8.854187817e-12
cimped = 376.730313461
cplanck = 6.626070040e-34 
cplanck_reduit = 1.054571800e-34 
cc = 299792458
cbolzman = 1.38064852e-23
cfaraday = 96485.33289
cStefBolz = 5.670367e-34 
cWien = 2.8977729e-3
cstructure = 7.2973525664

physics = [coulomb, ce, cG, cbolzman, 
    cfaraday, cstructure, cStefBolz, cWien, 
    cpermeability_magne, cpermeability_elec,
    cplanck, cplanck_reduit, cc, cimped]

constants = Array{Complex{Real}}([e, 1, 1im, pi, γ, φ, catalan])

operators = [*, +, /, -]
level = 0

function expr2string(expr::Expr)
    ope = expr.args[1]
    str_expr = "("

    for elem in expr.args[2:end-1]
        if typeof(elem) == Expr
            str_elem = expr2string(elem)
        else
            if imag(elem) == 0
                str_elem = String(string(real(elem)))
            elseif real(elem) == 0
                str_elem = String(string(imag(elem)))
            else
                str_elem = String(string(elem))
            end
        end
        str_expr = str_expr * str_elem * String(ope)
    end

    if typeof(expr.args[end]) == Expr
        last = expr2string(expr.args[end])
    else
        last = String(string(expr.args[end]))
    end

    return str_expr * last * ")"
end

function Base.hash(expr::Expr, h::UInt)
    base = hash((expr.head, expr.args[1]), h)
    if expr.args[1] in [:*, :+, :-]
        for expr in expr.args[2:end]
            base = base + hash(expr)
        end
    else
        base = base + hash(expr.args[2:end])
    end
    return base
end

function ==(a::Expr, b::Expr)     
    #return a === b || a.args[1] in [:+, :-, :*] && Set(a.args) == Set(b.args)
    if a.args[1] == :/ == b.args[1]
        return hash(a.args[2]) == hash(b.args[2]) && hash(a.args[3]) == hash(b.args[3])
    else 
        return hash(a) == hash(b)
    end
    
end

function reduce_add_mul(exp::Expr)
    for (i, elem) in enumerate(exp.args)
        if typeof(elem) == Expr && elem.args[1] == exp.args[1]
            deleteat!(exp.args, i)
            exp.args = vcat(exp.args, elem.args[2:end])
        end
    end
    return exp
end


function reduce_leaves(exp::Expr)
    for (i, elem) in enumerate(exp.args[2:end])
        if typeof(elem) == Expr
            if elem.args[1] in [:^, :/] && elem.args[3] == 1
                elem = elem.args[2]
            elseif elem.args[1] == :*
                args = filter(vals -> vals != 1, elem.args[2:end])
                if length(args) == 0
                    exp.args[i+1] = 1
                elseif length(args) == 1
                    exp.args[i+1] = args[1]
                else
                    exp.args[i+1] = Expr(:call, :*, args...) 
                end
            end
        end
    end

    if exp.args[1] == :*
        for (i, elem) in enumerate(exp.args[2:end])
            if typeof(elem) == Expr && elem.args[1] == :/
                to_add = exp.args[2:end]
                deleteat!(to_add, i)
                elem_left = vcat(elem.args[2], to_add)
                elem_right = elem.args[3]
                exp = Expr(:call, :/, Expr(:call, :*, elem_left...), elem_right)
                break
            end
        end
    end

    if exp.args[1] == :+
        for (i, elem) in enumerate(exp.args[2:end])
            if typeof(elem) == Expr && elem.args[1] == :/
                to_add = exp.args[2:end]
                deleteat!(to_add, i)
                elem_left = vcat(elem.args[2], to_add)
                elem_right = elem.args[3]
                exp = Expr(:call, :/, Expr(:call, :*, elem_left...), elem_right)
                break
            end
        end
    end

    return exp

end

function reduce_div2(expr::Expr)
    left = expr.args[2]
    right = expr.args[3]

    if typeof(left) == Expr && left.args[1] == :/
        left_left = left.args[2]
        left_right = left.args[3]
        expr = Expr(:call, :/, left_left, Expr(:call, :*,left_right, right))
    end

    if typeof(right) == Expr && right.args[1] == :/
        right_left = right.args[2]
        right_right = right.args[3]
        expr = Expr(:call, :/, Expr(:call, :*, left, right_right), right_left)
    end

    return expr

end

function reduce_div(exp::Expr)
    left = exp.args[2]
    right = exp.args[3]

    array_left = [left]
    array_right = [right]

    if typeof(right) == Expr
        if right.args[1] != :*
            return exp
        else        
            array_right = right.args[2:end]
        end
    end
    if typeof(left) == Expr
        if left.args[1] != :*
            return exp
        else        
            array_left = left.args[2:end]
        end
    end
    for (i, elem_left) in enumerate(array_left)
        for (j, elem_right) in enumerate(array_right)
            if length(array_left) == 0 || length(array_right) == 0
                break
            elseif elem_left == elem_right
                deleteat!(array_left, i)
                deleteat!(array_right, j)
            end
        end
    end
    
    if length(array_left) == 1
        new_left = array_left[1]
    elseif length(array_left) == 0
        new_left = 1 
    else
        new_left = Expr(:call, :*, array_left...)
    end

    if length(array_right) == 1
        new_right = array_right[1]
    elseif length(array_right) == 0
        new_right = 1
    else
        new_right = Expr(:call, :*, array_right...)
    end

    return Expr(:call, :/, new_left, new_right)
end

function rm_equi(itr)    
    seen = Vector{eltype(itr)}()
    for item in itr
        if all(x -> x != item, seen) 
            push!(seen, item)
        end
    end
    return seen
end

function isNullPotent(expr::Expr, null_potent)
    for elem in expr.args[2:end]
        if typeof(elem) == Expr && isNullPotent(elem, null_potent)
            return true
        end
    end
    return any(x -> x === expr, null_potent) 
end

function containSameCts(expr::Expr)
    all_args = expr.args[2:end]
    i = 1
    while i <= length(all_args)
        if typeof(all_args[i]) == Expr
            all_args = vcat(all_args, all_args[i].args[2:end])
            deleteat!(all_args, i)
        else
            i = i + 1
        end
    end
    return length(Set(all_args)) != length(all_args)
end

function reducing(exp::Expr)
    if exp.args[1] in [:+, :*]
        exp = reduce_add_mul(exp)
    elseif exp.args[1] == :/
        exp = reduce_div(reduce_div2(exp))
    end
    exp = reduce_leaves(exp)
    return exp
end

function applyFunct(cts)
    new_cts = Iterators.map(e -> [e, Expr(:call, :*, e, e), Expr(:call, :^, e, 0.5)], cts)
    return Iterators.vcat(new_cts...)
end

function applyOperator(left, operators, right)
    prod = Iterators.product(left, operators, right)
    return Iterators.map(x -> Expr(:call, Symbol(x[2]), x[1], x[3]), prod)
end

function oneExpr(cts, operators, next_level, null_potent, level)
    #println(typeof(cts), typeof(operators), typeof(next_level), typeof(null_potent), typeof(level))
    next_level2 = applyFunct(next_level)

    prod1 = applyOperator(cts, operators, next_level)
    prod2 = applyOperator(next_level, operators, cts)

    res = Iterators.vcat(prod1..., prod2...)
    res = Iterators.filter(expr -> !containSameCts(expr), res)
    res = Iterators.filter(expr -> !isNullPotent(expr, null_potent), res)
    res = collect(res)
    println(length(res))
    res = Iterators.map(x -> reducing(x), res)
    println(length(res))
    res = Iterators.unique(res)
    println(length(res))
    level = level + 1    
    #ee = Expr(:call, :*, e, Expr(:call, ))

    for ex in res
        r = eval(ex)
        if 0 == real(r) && abs(imag(r)) < 1.e-4
            push!(null_potent, ex)
        end
    end
    
    for exp in null_potent
        println(expr2string(exp), " ", eval(exp))
    end

    if level == 2
        return
    end

    oneExpr(cts, operators, res, null_potent, level)
end

the_good_one = Expr(:call, :*, pi + 0im, 1im)
the_good_one = Expr(:call, :^, e + 0im, the_good_one)
the_good_one = Expr(:call, :+, 1 + 0im, the_good_one)
#println(expr2string(the_good_one), " ", eval(the_good_one))


null_potent = Expr(:call, :*, 1im, 1im)
null_potent = Expr(:call, :^, 1im, null_potent)
test = Expr(:call, :^, 1im, null_potent)

reduce = Expr(:call, :*, 1im, pi + 0im)
reduce = Expr(:call, :/, pi + 0im, reduce)

reduce = Expr(:call, :/, pi + 0im, 1im)
reduce = Expr(:call, :/, pi + 0im, reduce)

reduce = Expr(:call, :/, 1im, pi + 0im)
reduce = Expr(:call, :/, pi + 0im, reduce)

#println(expr2string(reduce))
#reduce = reduce_div(reduce)
#println(expr2string(reduce))
#println(isNullPotent(test, [null_potent]))


null = Expr(:call, :/, 1, 1im)
null = Expr(:call, :+, 1im, null)



null2 = Expr(:call, :/, 1, 1im)
null2 = Expr(:call, :+, 1im, null2)

null = Expr(:call, :+, 1im, 1)


null2 = Expr(:call, :+, 1, 1im)




#println(expr2string(null), expr2string(null2), null==null2, null===null2)

constants = vcat(constants, math)
#names = vcat(names, math_names)

constants = vcat(constants, physics)
#names = vcat(names, physics_names)

#square(x) = x*x

#functions = [sqrt, square]

#labels = Dict(item[1] => item[2] for item in zip(constants, names))


test_diva = Expr(:call, :*, 1, 2)
test_divb = Expr(:call, :*, 2, 3, 6)
test_div = Expr(:call, :/, test_diva, test_divb)

test_div = reduce_leaves(test_div)

@time oneExpr(constants, operators, constants, Set(), level)

