include("../utils.jl")
include("../io.jl")
include("../graphs.jl")

@Logging.configure(level=INFO)

###
# loading and exporting datasets
###

# amazon_0601, web_google
dataset = "google"

#if dataset == "amazon"
	@info("loading Amazon_0601 graph")
	g = adjlist(UInt32, is_directed=true)
	load_adjacency_list_from_csv(UInt32, g, "../datasets/Amazon_0601/Amazon0601.txt")
	@info("# vertices:", length(vertices(g)))
	@info("# edges:", num_edges(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("# vertices:", length(vertices(core)))
	@info("# edges:", num_edges(core))

	@info("getting reverse graph")
	rcore = get_reverse_graph(core) 
	@info("# edges (rcore):", num_edges(rcore))

	write_mgs3_graph(core, "Amazon_0601_core")
	write_mgs4_graph(core, rcore, "Amazon_0601_core")
	serialize_to_jld(core, "core", "Amazon_0601_core")
	
	write_mgs3_graph(rcore, "Amazon_0601_rcore")
	write_mgs4_graph(rcore, core, "Amazon_0601_rcore")
	serialize_to_jld(rcore, "rcore", "Amazon_0601_rcore")
#elseif dataset == "google"
	@info("loading Web_Google graph")
	g = adjlist(UInt32, is_directed=true)
	load_adjacency_list_from_csv(UInt32, g, "../datasets/Web_Google/web-Google.txt")
	@info("# vertices:", length(vertices(g)))
	@info("# edges:", num_edges(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("# vertices:", length(vertices(core)))
	@info("# edges:", num_edges(core))

	@info("getting reverse graph")
	rcore = get_reverse_graph(core) 
	@info("# edges (rcore):", num_edges(rcore))

	write_mgs3_graph(core, "Web_Google_core")
	write_mgs4_graph(core, rcore, "Web_Google_core")
	serialize_to_jld(core, "core", "Web_Google_core")
	
	write_mgs3_graph(rcore, "Web_Google_rcore")
	write_mgs4_graph(rcore, core, "Web_Google_rcore")
	serialize_to_jld(rcore, "rcore", "Web_Google_rcore")
#end
