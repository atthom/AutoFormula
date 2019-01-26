using Base.MathConstants
using Profile
using BenchmarkTools
using Distributed
@everywhere using Reduce
using Base.Threads
using IterTools
using SharedArrays
#using AlgebraicNumbers

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

math = [cRS, cEB, csqrt2, cPlastique, cMills, cConway, cGK, 
    cApery, cKL, cBrun, cGauss, cLR, cET, clp, crobb,
    cpj, cC, cCD, comega, cHSM, cGKW, cbern, cmm,
    cKhintchine, cFeigenbaum, γ, φ, catalan]

math_names = ["cRS", "cEB", "csqrt2", "cPlastique", "cMills", "cConway", "cGK", 
    "cApery", "cKL", "cBrun", "cGauss", "cLR", "cET", "clp", "crobb",
    "cpj", "cC", "cCD", "comega", "cHSM", "cGKW", "cbern", "cmm",
    "cKhintchine", "cFeigenbaum", "γ", "φ", "catalan"]

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

constants = Array{Real}([e, 1, pi])
cts_names = Array{String}(["e", "1", "pi", "1im"])
#i removed

operators = ["*", "+", "/", "-"]
level = 0

function applyOperator(left, operators, right)
    prod = Iterators.product(left, operators, right)
    return map(x -> Expr(:call, Symbol(x[2]), Symbol(x[1]), Symbol(x[3])), prod)
end

function expr2file(expr)
    #return String(string(expr)) * "\n"
    return string(expr) * "\n"
end

function save_to_file(array, file_name)
    if isempty(array)
        println("Not saving, null_potent is empty...")
    else
        println("Saving expressions...")
        str_array = Iterators.mapreduce(expr2file, *, array)
        write(open(file_name * ".txt", "w"), str_array)
    end
end

function keep_null_potent(res)
    null_potent = filter(item -> typeof(item[2])==Int && item[2] == 0, res)
    println("Null potent length:", length(null_potent))
    save_to_file(map(item -> item[1], null_potent), "null_potent")
end

function oneExpr(cts, operators, next_level, level)
    level = level + 1
    #new_cts = Iterators.map(e -> [e, Expr(:call, :*, e, e), Expr(:call, :^, e, 0.5)], cts)
    #new_cts = Iterators.vcat(new_cts...)

    prod1 = applyOperator(cts, operators, next_level) # 5.160 ms 448.203 ms
    prod2 = applyOperator(next_level, operators, cts) # 5.198 ms 402.022 ms

    res = vcat(prod1..., prod2...) # 31.118 μs 9.209 ms
    println(length(res))
    
    #@distributed res = map(exp -> (exp, rcall(exp)), res)
    
    #@everywhere new_res = []
    #@sync @distributed res for exp in res
    #    push!(new_res, [exp, rcall(exp)])
    #end
    #new_res = SharedArray{Expr}(length(res))
    #@sync @distributed for exp in res
    #    push!(new_res, rcall(exp))
    #end
    res = map(exp -> (exp, rcall(exp)), res)
    
    #@everywhere new_res = []
    #@sync @distributed for i in 1:length(res)
    #    push!(new_res, (res[i], rcall(res[i])))
    #end
    #res = new_res
    #@time res = pmap(exp -> (exp, rcall(exp)), res, distributed=false)
    @time keep_null_potent(res)
    @time res = map(item -> item[2], res)
    @time res = unique(filter(exp -> typeof(exp) == Expr, res))
    println("result set length:", length(res))

    if level == 3
        return
    end

    @time oneExpr(cts, operators, res, level)
end

constants = vcat(constants, math)
cts_names = vcat(cts_names, math_names)

#constants = vcat(constants, physics)
#cts_names = vcat(cts_names, physics_names)
labels = Dict(item[1] => item[2] for item in zip(constants, cts_names))

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    oneExpr(cts_names, operators, cts_names, Set(), level)
    return 0
end

@time oneExpr(cts_names, operators, cts_names, level)
