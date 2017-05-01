include("../utils.jl")
include("../io.jl")
include("../graphs.jl")
include("../PR.jl")
include("../RW.jl")

@Logging.configure(level=INFO)

###
# shortest distance
###

@info("############ Testing shortest paths functions")

core = adjlist(UInt32, is_directed=true)
load_mgs3_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")
#load_mgs4_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgz")

@info("getting rcore")
rcore = get_reverse_graph(core)

@info("transforming core from adj to inc list")
core2 = inclist(vertices(core), is_directed=true)

@info("adding edges to inc list")
for u in vertices(core)
	for v in out_neighbors(u,core)
		add_edge!(core2,u,v)
	end
end
#serialize_to_jld(core2, "core", "Arxiv_HEP-PH_core_inclist")

@info("computing Djikstra from vertex 1")
dists = ones(num_edges(core2))
s = convert(UInt32,1)
t = convert(UInt32,100)
@time r = dijkstra_shortest_paths(core2, dists, s)

@info("mean distance from vertex 1: ", mean(r.dists))
paths = enumerate_paths(vertices(core2), r.parent_indices, [10,100])
@info("shortest paths from vertex 1 to vertices 10 and 100: ", convert(Array{Array{Int64,1},1}, paths))

@info("computing Bellman-Ford from vertex 1")
@time r = bellman_ford_shortest_paths(core2, dists, UInt32[1])

@info("mean distance from vertex 1: ", mean(r.dists))
paths = enumerate_paths(vertices(core2), r.parent_indices, [10,100])
@info("shortest paths from vertex 1 to vertices 10 and 100: ", convert(Array{Array{Int64,1},1}, paths))

@info("computing A* from node 1 to vertex 100")
@time r = shortest_path(core2, dists, s, t)
@info("shortest path from vertex 1 to vertex 100: ", r)

@info("computing Pagerank of core and rcore")
eps = 1e-4
@time pr_core = PR(core, rcore, epsilon=eps)
@time pr_rcore = PR(rcore, core, epsilon=eps)

@info("computing personalized Pageranks for nodes 1 and 100")
@time pr_s = personalized_PR(s, core, rcore, epsilon=eps)
@time pr_t = personalized_PR(t, rcore, core, epsilon=eps)

@info("computing diffusion fingerprints")
df = pr_s .* pr_t
pr_boost = pr_core .* pr_rcore
df_boost = df ./ pr_boost

# max length
ml = log(length(vertices(core)))
# current vertex
cv = s
# current path
sp = UInt32[]

@info("finding shortest path (greedy approach)")

# greedy approach
while length(sp) < ml
	@debug("adding vertex: ", cv)
	push!(sp,cv)
	nei = out_neighbors(cv,core)
	nnei = setdiff(nei,sp)
	if length(nnei) == 0
		@info("--- exploration reached a dead end")
		@info("--- explored path: ", sp)
		break
	end
	cv = nnei[indmax(df[nnei])]
	#cv = nei[indmax(df_boost[nnei])]

	# we found a path
	if cv == t
		@debug("path found: ", sp)
		push!(sp,t)
		break
	end
end

@info("shortest path: ", convert(Array{Int64,1}, sp))
@info("finding shortest path (probabilistic approach)")

# array of paths
sps = Array{Array{UInt32,1},1}()
# current path
sp = UInt32[]
cv = s

# search in a probabilistic way shortest paths between s and t
c = 0
max_iter = 1e3
while c < max_iter
	while length(sp) < ml 
		@debug("adding vertex: ", cv)
		push!(sp,cv)
		nei = out_neighbors(cv,core)
		nnei = setdiff(nei,sp)
		if length(nnei) == 0
			@debug("--- exploration reached a dead end")
			@debug("--- explored path: ", sp)
			sp = UInt32[]
			cv = s
			break
		end
		pos = get_flying_index(df[nnei] / sum(df[nnei]))
		cv = nnei[pos]

		# we found a path
		if cv == t
			@debug("path found: ", sp)
			push!(sp,t)
			if !(sp in sps)
				push!(sps,sp)
			end
			sp = UInt32[]
			cv = s
			break
		end
	end
	c += 1
end

@info("shortest paths: ", convert(Array{Array{Int64,1},1}, sps))
