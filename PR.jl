#
# JCNL: Julia Complex Networks Library
# Copyright (C) 2016-2017  Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

include("RW.jl")

using Graphs, DataStructures, Logging, Distances

# simple implementation of Pagerank algorithm
function PR{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}; init_pr::Array{Float64,1}=Float64[], damping::Float64=0.85, epsilon::Float64=1e-4)
	vs = vertices(g)
	n = length(vs)
	@info("computing Pagerank (size of graph $n)")
	# initialize pagerank vector
	if length(init_pr) == n
		pr = init_pr
	else
		pr = Float64[1/n for k in vs]
	end
	pr2 = copy(pr)
	while true
		for v in vs
			nv = 0.
			# get v children in the reverse graph
			in_nei = out_neighbors(v,rg)
			if length(in_nei) > 0
				for p in in_nei
					nv +=  pr[p]/length(out_neighbors(p,g))
				end
			end
			pr2[v] = (1-damping)/n+damping*nv
		end
		d = chebyshev(pr,pr2)
		pr = copy(pr2)
		d <= epsilon && break
	end
	return pr
end

# simple implementation of Pagerank algorithm
#
# specifically designed for large graphs
function PR{T<:Unsigned}(rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, out_degrees::Array{T,1}; init_pr::Array{Float64,1}=Float64[], damping::Float64=0.85, epsilon::Float64=1e-4, save_pr::Bool=False)
	vs = vertices(rg)
	n = length(vs)
	@info("computing Pagerank (size of graph $n)")
	# initialize pagerank vector
	# initialize pagerank vector
	if length(init_pr) == n
		pr = init_pr
	else
		pr = Float64[1/n for k in vs]
	end
	pr2 = copy(pr)	
	# iteration number
	ic = 1
	# progression
	thp = ceil(Int,n/100)
	while true
		for v in vs
			nv = 0.
			# get v children in the reverse graph
			in_nei = out_neighbors(v,rg)
			if length(in_nei) > 0
				for p in in_nei
					nv +=  pr[p]/out_degrees[p]
				end
			end
			pr2[v] = (1-damping)/n+damping*nv
			if v % thp == 0
				# log progression
				@info("iteration $ic: ", (v/n)*100, " %")
			end
		end
		d = chebyshev(pr,pr2)
		pr = copy(pr2)
		@info("distance(t,t+1) = $d")
		if save_pr
			serialize_to_file(pr, "pr-iter-$ic.jld")
		end
		ic = ic + 1
		d <= epsilon && break
	end
	return pr
end

# simple implementation of personalized Pagerank with a single source vertex
function personalized_PR{T<:Unsigned}(src::T, g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}; damping::Float64=0.85, epsilon::Float64=1e-4)
	vs = vertices(g)
	n = length(vs)
	@info("computing personalized Pagerank (size of graph $n, source $src)")
	# initialize pagerank vector
	pr = [0. for k in vs]
	pr[src] = 1.
	pr2 = copy(pr)
	while true
		for v in vs
			nv = 0.
			# get v children in the reverse graph
			in_nei = out_neighbors(v,rg)
			if length(in_nei) > 0
				for p in in_nei
					nv +=  pr[p]/length(out_neighbors(p,g))
				end
			end
			if v == src
				pr2[v] = (1-damping)+damping*nv
			else
				pr2[v] = damping*nv
			end
		end
		d = chebyshev(pr,pr2)
		pr = copy(pr2)
		d <= epsilon && break
	end
	return pr
end


# Monte-Carlo Pagerank algorithm
#
# MC Pagerank with cyclic start of complete path stopping at sink nodes
# http://www-sop.inria.fr/members/Konstantin.Avratchenkov/pubs/mc.pdf
function MC_PR{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, n_cycles::Int; damping::Float64=0.85, epsilon::Float64=1e-4)
	vs = vertices(g)
	n = length(vs)
	vv = zeros(Float64, n)
	@info("computing Monte-Carlo Pagerank (size of graph $n)")
	for i in 1:n_cycles
		for v in vs
			vv += RW_aggregated(g, (1-damping), v)
		end
	end
	return vv/sum(vv)
end
