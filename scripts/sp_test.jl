include("../utils.jl")
include("../io.jl")
include("../graphs.jl")

@Logging.configure(level=INFO)

###
# shortest distance
###

@info("############ Testing shortest paths functions")

core = adjlist(UInt32, is_directed=true)
load_mgs3_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")
#load_mgs4_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgz")

@info("transforming core from adj to inc list")
core2 = inclist(vertices(core), is_directed=true)

@info("adding edges to inc list")
for u in vertices(core)
	for v in out_neighbors(u,core)
		add_edge!(core2,u,v)
	end
end
#serialize_to_jld(core2, "core", "Arxiv_HEP-PH_core_inclist")

dists = ones(num_edges(core2))

@info("computing Djikstra from vertex 1")
@time r = dijkstra_shortest_paths(core2, dists, convert(UInt32,1))

@info("mean distance from vertex 1: ", mean(r.dists))
paths = enumerate_paths(vertices(core2), r.parent_indices, [10,100])
@info("paths from vertex 1 to vertices 10 and 100: ", convert(Array{Array{Int64,1},1}, paths))

@info("computing Bellman-Ford from vertex 1")
@time r = bellman_ford_shortest_paths(core2, dists, UInt32[1])

@info("mean distance from vertex 1: ", mean(r.dists))
paths = enumerate_paths(vertices(core2), r.parent_indices, [10,100])
@info("paths from vertex 1 to vertices 10 and 100: ", convert(Array{Array{Int64,1},1}, paths))

@info("computing A* from node 1 to vertex 100")

@time r = shortest_path(core2, dists, convert(UInt32,1), convert(UInt32,100))
@info("path: ", r)
