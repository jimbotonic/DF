include("../graphs.jl")
include("../io.jl")
include("../PR.jl")

@Logging.configure(level=DEBUG)

filename = ARGS[1]

# load graph in MGSv3 format
g = adjlist(UInt32, is_directed=true)
load_mgs3_graph(g, filename)
rg = get_reverse_graph(g)

tic()
pr = PR(g, rg)
toc()

tic()
pr_mc = MC_PR(g, 200)
toc()

println(pr[1:20])
println(pr_mc[1:20])
println(chebyshev(pr, pr_mc))
