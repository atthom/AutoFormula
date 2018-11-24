using Combinatorics

@inline function dist(p1::Tuple{Int64, Int64}, p2::Tuple{Int64, Int64})::Int64
    return (p2[1] - p1[1])^2 + (p2[2] - p1[2])^2;
end

@inline function dist(cities::Array{Tuple{Int64, Int64}})::Int64
    cities2 = copy(cities[2:end]);
    push!(cities2, cities[1]);
    return sum(map(p -> dist(p[1], p[2]), zip(cities, cities2)));
end

function traveling(cities::Array{Tuple{Int64, Int64}})
    all_permutations = permutations(cities[2:end])
    println(length(all_permutations), " permutations")

    all_permutations = map(perm -> pushfirst!(perm, cities[1]), all_permutations);

    all_distances = map(perm -> (perm, dist(perm)), all_permutations);
    result = sort(all_distances, by = x -> x[2]);
    println(result[1])
end

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    cities = [(1, 2), (14, 27), (15, 19), (4, 9), (0, 50), (5, 13), (13, 22), (13, 22), (25, 6), (30, 14), (7, 7)];
    #cities = [(1, 2), (14, 27), (4, 9), (0, 50), (5, 13)];
    @time traveling(cities);
    return 0
end
#cities = [(1, 2), (14, 27), (15, 19), (4, 9), (0, 50), (5, 13), (13, 22), (13, 22), (25, 6), (30, 14), (7, 7)];
#@time traveling(cities);

julia_main(Vector{String}([]))