include("../io.jl")
include("../graphs.jl")

@Logging.configure(level=DEBUG)

g = adjlist(UInt32, is_directed=true)

filename = ARGS[1]

# load CSV adjacency list
load_graph_from_pajek(UInt32, g, filename)

nvs,nes,dens = get_basic_stats(g)
println(nvs, " ", nes," ", dens)

# save graph in MGSv3 format
write_mgs3_graph(g, "EAT")

# load graph in MGSv3 format
g2 = adjlist(UInt32, is_directed=true)
load_mgs3_graph(g2, "EAT.mgs")

nvs,nes,dens = get_basic_stats(g2)
println(nvs, " ", nes," ", dens)

vs = vertices(g2)
for v in vertices(g2)[1:1000]
	println(out_neighbors(v,g2))
end
