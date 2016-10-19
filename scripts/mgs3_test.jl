include("../io.jl")
include("../graphs.jl")

@Logging.configure(level=DEBUG)

g = adjlist(UInt32, is_directed=true)

filename = ARGS[1]

# load CSV adjacency list
load_adjacency_list_from_csv(UInt32, g, filename, '\t')

nvs,nes,dens = get_basic_stats(g)

# save graph in MGSv3 format
write_mgs3_graph(g, "Arxiv_HEP-PH")

# load graph in MGSv3 format
g2 = adjlist(UInt32, is_directed=true)
load_mgs3_graph(g2, "Arxiv_HEP-PH.mgs")

println(length(vertices(g)))
println(out_neighbors(vertices(g)[1],g))
println(length(vertices(g2)))
println(out_neighbors(vertices(g2)[1],g2))
#println(length(edges(g2)))
