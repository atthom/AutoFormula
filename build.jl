using PackageCompiler

#build_executable("TravelingSalesman.jl")

dd = normpath(Base.find_package("PackageCompiler"), "..", "..")
print(dd)