#
# JCNL: Julia Complex Networks Library
# Copyright (C) 2016  Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

using Graphs, DataStructures, Logging

# naïve implementation of Pagerank algorithm
function my_pagerank{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},pr::Array{Float64,1},damping::Float64=0.85,epsilon::Float64=1e-4)
	vs = vertices(g)
	n = length(vs)
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
		d = dist(pr,pr2)
		pr = copy(pr2)
		d <= epsilon && break
	end
	return pr
end

function my_pagerank2{T<:Unsigned}(rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},out_degrees::Array{T,1},damping::Float64=0.85,epsilon::Float64=1e-4,save_pr::Bool=False,init_pr=None)
	vs = vertices(rg)
	n = length(vs)
	@debug("computing pagerank (size of graph $n)")
	# initialize pagerank vector
	if init_pr == None
		pr = [float64(1/n) for k in vs]
	else
		pr = init_pr
	end
	pr2 = copy(pr)
	ic = 1
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
			if v % 100000 == 0
				@debug("iteration $ic: ", (v/n)*100, " %")
			end
		end
		d = dist(pr,pr2)
		pr = copy(pr2)
		@debug("distance(t,t+1) = $d")
		if save_pr
			serialize_to_file(pr, "pr-iter-$ic.jld")
		end
		ic = ic + 1
		d <= epsilon && break
	end
	return pr
end

# naïve implementation of personalized Pagerank with a single source vertex
function my_personalized_pagerank{T<:Unsigned}(src::T,g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},damping::Float64=0.85,epsilon::Float64=1e-4)
	vs = vertices(g)
	n = length(vs)
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
		d = dist(pr,pr2)
		pr = copy(pr2)
		d <= epsilon && break
	end
	return pr
end
