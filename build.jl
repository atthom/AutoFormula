#using PackageCompiler

#build_executable("TravelingSalesman.jl")

#dd = normpath(Base.find_package("PackageCompiler"), "..", "..")
#print(dd)

#build_executable("GenerateFormula3.jl")

path = normpath(Base.find_package("PackageCompiler"), "..", "..")
command = "julia " * path * "juliac.jl -vRe ./GenerateFormula3.jl"
#print(command)
run(`$command`)
