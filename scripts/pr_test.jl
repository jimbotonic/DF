include("../graphs.jl")
include("../io.jl")
include("../pr.jl")

# load graph in MGSv3 format
core = adjlist(UInt32, is_directed=true)
load_mgs3_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")
#load_mgs4_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgz")

@info("getting rcore")
rcore = get_reverse_graph(core)

# constants
n = num_vertices(core)
damping = .85
eps = 1e-8

@info("computing Pagerank of core and rcore")
@time pr_core = PR(core, rcore, epsilon=eps)
#@time pr_rcore = PR(rcore, core, epsilon=eps)

@info("computing personalized Pageranks for nodes 1 and 100")
# source and target node
s = convert(UInt32,1)
t = convert(UInt32,100)
@time pr_s = PPR(s, core, rcore, epsilon=eps)
#@time pr_t = PPR(t, rcore, core, epsilon=eps)

## compute Monte Carlo PR
#niter = 400
#@info("compute Pagerank (Monte Carlo)")
#@time pr_mc = MC_PR(core, niter)
#
#@info("pr_core <-> pr_mc: ", chebyshev(pr_core, pr_mc))

## matrix computation
@info("getting P matrix")
P = get_sparse_P_matrix(core)

@info("compute Pagerank (power iteration)")
ppr = zeros(Float64,n)
ppr[s] = 1.
@time pr_pi = PR(P, epsilon=eps)

@info("computing personalized Pageranks for node 1 (power iteration)")
@time pr_pi_s = PR(P, ppr=ppr, epsilon=eps)

@info("pr_core <-> pr_pi: ", chebyshev(pr_core, pr_pi))
@info("pr_s <-> pr_pi_s: ", chebyshev(pr_s, pr_pi_s))
