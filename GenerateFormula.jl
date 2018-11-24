
using Base.MathConstants

import Base.==

constants = Array{Complex{Real}}([e, pi, 1, 1im])
names = Array{String}(["e", "π", "1", "i"])
labels = Dict(item[1] => item[2] for item in zip(constants, names))

functions = [^, *, +]
level = 0

function eq_modulo_ordering!(xs, ys, comm)  # note !, mutates xs and ys
    while !isempty(xs)
        x = pop!(xs)
        i = findfirst(b -> expr_equiv(x, b, comm), ys)
        i === nothing && return false
        deleteat!(ys, i)
    end
    isempty(ys)
end
eq_modulo_ordering(xs, ys, comm) = eq_modulo_ordering!(copy(xs), copy(ys), comm)

function expr_equiv(a::Expr, b::Expr, comm)
    a.head === b.head || return false
    a.head === :call || return a == b
    a.args[1] in comm || return all(expr_equiv.(a.args, b.args, Ref(comm)))

    eq_modulo_ordering(a.args, b.args, comm)
end
expr_equiv(a, b, comm) = a == b
expr_equiv(a, b) = expr_equiv(a, b, [:+, :*])

function ==(a::Expr, b::Expr) 
    return expr_equiv(a, b)
end



function expr2string(expr)
    ope = expr.args[1]
    left = expr.args[2]
    right = expr.args[3]
    if typeof(left) == Expr
        left = expr2string(left)
    else
        left = labels[left]
    end
    if typeof(right) == Expr
        right = expr2string(right)
    else
        right = labels[right]
    end
    
    return "(" * left * String(ope) * right * ")"
end

function remove_equivalences(itr)
    new_itr = Vector{Expr}([])
    for expr in itr
        if all(item -> item != expr, new_itr)
            push!(new_itr, expr)
        end
    end
    #new_itr = Set{Expr}(itr)
    #println(length(itr), " ", length(new_itr))
    return new_itr
end


function rm_equi!(itr)    
    seen = Set{eltype(itr)}()
    for (i, item) in enumerate(itr)
        if any(x -> x == item, seen) 
            deleteat!(itr, i)
        else
            push!(seen, item)
        end
    end
    return itr
end

function isNullPotent(expr, null_potent)
    if in(expr, null_potent)
        return true
    elseif expr.args[2] == Expr && isNullPotent(expr.args[2], null_potent)
        return true
    elseif expr.args[3] == Expr && isNullPotent(expr.args[3], null_potent)
        return true
    end
    return false
end



function oneExpr(cts, functions, next_level, null_potent, level)
    #next_level = vcat(next_level, cts)
    prod = Iterators.product(cts, functions, next_level)
    res = Iterators.map(x -> Expr(:call, Symbol(x[2]), x[1], x[3]), prod)
    res = Iterators.vcat(res...)
    println(size(res))
    rm_equi!(res)
    println(size(res))
    res = Iterators.filter(expr -> !isNullPotent(expr, null_potent), res)
    res = Iterators.filter(expr -> expr.args[1]==Symbol(+), res)
    res = collect(res)
    println(size(res))

    level = level + 1

    for ex in res
        r = eval(ex)
        if 0 == real(r) && imag(r) < 1.e-4
            push!(null_potent, ex)
        end
    end

    for exp in null_potent
        println(expr2string(exp))
    end

    for exp in res
        println(expr2string(exp))
    end

    if level == 2
        return
    end

    oneExpr(cts, functions, res, null_potent, level)
end

the_good_one = Expr(:call, :*, pi + 0im, 1im)
the_good_one = Expr(:call, :^, e + 0im, the_good_one)
the_good_one = Expr(:call, :+, 1 + 0im, the_good_one)
#println(expr2string(the_good_one), " ", eval(the_good_one))

@time oneExpr(constants, functions, constants, [], level)

