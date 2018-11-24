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

math_names = ["cRS", "cEB", "cGolden", "csqrt2", "cPlastique", "cMills", "cConway", "cGK", 
    "cApery", "cKL", "cK", "cBrun", "cGauss", "cLR", "cET", "clp", "crobb",
    "cpj", "cC", "cCD", "cEM", "comega", "cHSM", "cGKW", "cbern", "cmm",
    "cKhintchine", "cFeigenbaum"]

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

physics_names = ["coulomb", "ce", "cG", "cbolzman", 
    "cfaraday", "cstructure", "cStefBolz", "cWien", 
    "cpermeability_magne", "cpermeability_elec",
    "cplanck", "cplanck_reduit", "cc","cimped"]

constants = Array{Complex{Real}}([e, 1, pi, 1im, γ, φ, catalan])
names = Array{String}(["e", "1", "π", "i", "γ", "φ", "catalan"])

operators = [*, +, /]
level = 0

function ==(a::Expr, b::Expr)     
    return a === b || a.args[1] in [:+, :*, :/] && Set(a.args) == Set(b.args)
end

function reduce_add_mul(exp)
    left = exp.args[2]
    if typeof(left) == Expr && left.args[1] in [:+, :*]
        exp.args = vcat(exp.args, left.args[2:end])
        deleteat(exp.args, 2)
    end
    right = exp.args[3]
    if typeof(right) == Expr && right.args[1] in [:+, :*]
        exp.args = vcat(exp.args, right.args[2:end])
        deleteat!(exp.args, 3)
    end
end

function reduce_power(exp)
    left = exp.args[2]
    if typeof(left) == Expr && left.args[1] == :^ && left.args[3] == 1
        exp.args[2] = left.args[2]
    end
    right = exp.args[3]
    if typeof(right) == Expr && right.args[1] == :^ && right.args[3] == 1
        exp.args[3] = right.args[2]
    end
    return exp
end

function reduce_div(exp)
    left = exp.args[2]
    right = exp.args[3]
    
    if exp.args[1] == :/ && typeof(right) == Expr

        if right.args[1] == :/ && left == right.args[2]
            return Expr(:call, :*, 1,  right.args[3])
        elseif left == right.args[3]
            if right.args[1] == :*
                return Expr(:call, :/, 1,  right.args[2])
            elseif right.args[1] == :/
                return Expr(:call, :*, left, 1)
            end
        elseif right.args[1] == :* && left == right.args[2]
            return Expr(:call, :/, 1, right.args[3])
        end
    end

    return exp
end

function reducing(exp::Expr)
    if exp.args[1] in [:+, :*]
        reduce_add_mul(exp)
    end
    exp = reduce_power(exp)
    exp = reduce_div(exp)
    return exp
end


function expr2string(expr)
    ope = expr.args[1]
    str_expr = "("

    for elem in expr.args[2:end-1]
        if typeof(elem) == Expr
            str_elem = expr2string(elem)
        else
            str_elem = labels[elem]
        end
        str_expr = str_expr * str_elem * String(ope)
    end

    if typeof(expr.args[end]) == Expr
        last = expr2string(expr.args[end])
    else
        last = labels[expr.args[end]]
    end

    return str_expr * last * ")"
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


function isNullPotent(expr, null_potent)
    if typeof(expr.args[2]) == Expr && isNullPotent(expr.args[2], null_potent)
        return true
    elseif typeof(expr.args[3]) == Expr && isNullPotent(expr.args[3], null_potent)
        return true
    end
    return any(x -> x == expr, null_potent) 
end

function oneExpr(cts, operators, next_level, null_potent, level)
    prod = Iterators.product(cts, operators, next_level)
    res = Iterators.map(x -> Expr(:call, Symbol(x[2]), x[1], x[3]), prod)
    res = Iterators.vcat(res...)

    prod2 = Iterators.product(next_level, operators, cts)
    res2 = Iterators.map(x -> Expr(:call, Symbol(x[2]), x[1], x[3]), prod2)
        
    res = Iterators.vcat(res..., res2...)
    #res = Iterators.map(x -> reducing(x), res)
    println(length(res))
    res = Iterators.unique(res)
    println(length(res))
    #res = rm_equi(res)
    println(length(res))
    res = Iterators.filter(expr -> !isNullPotent(expr, null_potent), res)
    res = collect(res)
    println(length(res))

    level = level + 1

    for ex in res
        r = eval(ex)
        if 0 == real(r) && abs(imag(r)) < 1.e-4
            push!(null_potent, ex)
        end
    end

    for exp in null_potent
        println(expr2string(exp), " ", eval(exp))
    end

    if level == 3
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

#println(expr2string(null), expr2string(null2), null==null2, null===null2)

for (i, elem) in enumerate(math)
    #push!(constants, elem)
    #push!(names, math_names[i])
end

for (i, elem) in enumerate(physics)
    #push!(constants, elem)
    #push!(names, physics_names[i])
end

#constantes_squared = []
#constantes_squared_names = []
#for (i, constant) in enumerate(constants)
#    push!(constantes_squared_names, names[i] * "²")
#    push!(constantes_squared, constant*constant)
#end

#for (i, constant) in enumerate(constantes_squared)
#    push!(constants, constant)
#    push!(names, constantes_squared_names[i])
#end

labels = Dict(item[1] => item[2] for item in zip(constants, names))

@time oneExpr(constants, operators, constants, Set(), level)

