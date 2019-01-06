
the_good_one = Expr(:call, :*, pi + 0im, 1im)
the_good_one = Expr(:call, :^, e + 0im, the_good_one)
the_good_one = Expr(:call, :+, 1 + 0im, the_good_one)
#println(expr2string(the_good_one), " ", eval(the_good_one))

null_potent = Expr(:call, :*, 1im, 1im)
null_potent = Expr(:call, :^, 1im, null_potent)
test = Expr(:call, :^, 1im, null_potent)

red = Expr(:call, :*, 1im, pi + 0im)
red = Expr(:call, :/, pi + 0im, red)

red = Expr(:call, :/, pi + 0im, 1im)
red = Expr(:call, :/, pi + 0im, red)

red = Expr(:call, :/, 1im, pi + 0im)
red = Expr(:call, :/, pi + 0im, red)

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

#square(x) = x*x

#functions = [sqrt, square]


test_diva = Expr(:call, :*, 1, 2)
test_divb = Expr(:call, :*, 2, 3, 6)
test_div = Expr(:call, :/, test_diva, test_divb)
#test_div = reduce_leaves(test_div)
#constants = map(ct -> AlgebraicNumber(ct), cts)
#println(expr2string(null), expr2string(null2), null==null2, null===null2)