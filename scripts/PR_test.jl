include("../graphs.jl")
include("../io.jl")
include("../PR.jl")

@Logging.configure(level=DEBUG)

filename = ARGS[1]

# load graph in MGSv3 format
g = adjlist(UInt32, is_directed=true)
load_mgs3_graph(g, filename)

tic()
core,oni,noi = get_core(g)
toc()


tic()
rcore = get_reverse_graph(core)
toc()

write_mgs3_graph(core, "Arxiv_HEP-PH_core")

tic()
pr = PR(core, rcore, epsilon=1e-8)
toc()

serialize_to_jld(pr, "PR", "Arxiv_HEP-PH_PR")

tic()
pr_mc = MC_PR(core, 400)
toc()

println(pr[1:20])
println(pr_mc[1:20])
println(chebyshev(pr, pr_mc))
