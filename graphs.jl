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

###
# Basic stats for directed graphs
###

# get # vertices, # of edges and density
# NB: for computing density, we assume a directed graph
function get_basic_stats(g)
	nvs = length(vertices(g))
	nes = num_edges(g)
	density = nes/(nvs*(nvs-1))
	return nvs,nes,density
end

# display basic stats
function display_basic_stats(g,rg)
	nvs,nes,dens = get_basic_stats(g)
	nvs2,nes2,dens2 = get_basic_stats(rg)

	avg_od,max_od,sinks = get_out_degree_stats(g)
	avg_id,max_id,sources = get_out_degree_stats(rg)

	@info("# vertices: $nvs")
	@info("# edges: $nes")
	@info("density: $dens")
	@info()
	@info("avg out-degree: $avg_od")
	@info("max out-degree: $max_od")
	nsinks = length(sinks)
	@info("# sinks: $nsinks")
	@info()
	@info("avg in-degree: $avg_id")
	@info("max in-degree: $max_id")
	nsources = length(sources)
	@info("# sources: $nsources")
end

# get out degree stats
#
# retuns avg out-degree, max out-degree, array of sinks vertices
function get_out_degree_stats{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	sum = 0
	sinks = T[]
	max_degree = 0
	vs = vertices(g)
	for v in vs
		children = out_neighbors(v,g)
		od = length(children)
		sum += od
		if od == 0
			push!(sinks,v)
		end
		if od > max_degree
			max_degree = od
		end
	end
	return sum/length(vs),max_degree,sinks
end

# get array of sink vertices in the specified graph
function get_sinks{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	sinks = T[]
	vs = vertices(g)
	for v in vs
		children = out_neighbors(v,g)
		od = length(children)
		if od == 0
			push!(sinks,v)
		end
	end
	return sinks
end

# get array of source vertices in the specified graph
function get_sources{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	achildren = T[]
	vs = vertices(g)
	for v in vs
		children = out_neighbors(v,g)
		achildren = union(achildren,children)
	end
	return setdiff(vs,achildren)
end

###
# Subgraphs && SCCs
###

# get the subgraph of g induced by the set of vertices sids
function subgraph{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},sids::Array{T,1})
	ng = adjlist(T, is_directed=true)
	# nvs should be sorted in acending order
	nvs = sort(sids)

	# add vertices
	for i in 1:length(nvs)
        	add_vertex!(ng,convert(T,i))
	end

	# old id -> new id
	oni = Dict{T,T}()
	noi = Dict{T,T}()

	counter = convert(T,1)
	for v in nvs
		oni[v] = counter
		noi[counter] = v
		counter += convert(T,1)
	end

	# add edges
	for v in nvs
		children = out_neighbors(v,g)
		source = oni[v]
		for c in children
			#if index_sorted(nvs,c) != 0
			if haskey(oni,c)
				target = oni[c]
				add_edge!(ng,source,target)
			end
		end
	end
	return ng,oni,noi
end

# get the subgraph of g induced by the set of vertices sids
# write the subgraph in MGS format v2
function subgraph_streamed{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, sids::Array{T,1}, name::String)
	ng = adjlist(T, is_directed=true)
	# nvs should be sorted in acending order
	nvs = sort(sids)

	f1 = open("$name.index", "w")
	f2 = open("$name.data", "w")

	# old id -> new id
	oni = Dict{T,T}()

	counter = convert(T,1)
	for v in nvs
		oni[v] = counter
		counter += convert(T,1)
	end

	pos = convert(T,1)
	for v in nvs
		children = out_neighbors(v,g)
		bytes = reinterpret(Uint8, [pos])
                write(f1, reverse(bytes))
		for c in children
			if haskey(oni,c)
				target = oni[c]
				bytes = reinterpret(Uint8, [target])
                		write(f2, reverse(bytes))
				# pos += convert(Uint32, 1) -> 8 bytes!?
				pos = convert(T,pos+1)
			end
		end
	end

	@debug("# vertices: ", length(nvs))
	@debug("# edges: ", pos)

	close(f1)
	close(f2)

	return oni
end

# Tarjan algorithm (recursive version)
#
# NB: successfully tested with FA core
# NB: the recursive calls may create a stack overflow error
function tarjan{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	sccs = Array(Array{T,1},0)
	n = length(vertices(g))
	indices = zeros(T,n)
	lowlinks = zeros(T,n)
	S = T[]
	index = 1

	# recursive function
	function visit(v)
		indices[v] = index
		lowlinks[v] = index
		index += 1
		push!(S,v)
		children = out_neighbors(v,g)
		for w in children
			# w was not visited yet
			if indices[w] == 0
				# recursive call
				visit(w)
				# back propagate the lowlink
				lowlinks[v] = min(lowlinks[v], lowlinks[w])
			# w is in S and hence in the current SCC
			elseif findfirst(S,w) != 0
				lowlinks[v] = min(lowlinks[v], indices[w])
			end
		end

		if lowlinks[v] == indices[v]
			component = T[]
			while true
				w = pop!(S)
				push!(component,w)
				w == v && break
			end
			push!(sccs,component)
		end
	end

	for v in vertices(g)
		indices[v] == 0 && visit(v)
	end

	return sccs
end

# Pearce version of Tarjan algorithm
#
# NB: successfully tested with FA core
# NB: the recursive calls may create a stack overflow error
function pearce{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	n = length(vertices(g))
	rindex = zeros(T,n)
	S = T[]
	index = 1
	c = n-1

	# recursive function
	function visit(v)
		root = true
		rindex[v] = index
		index += 1
		children = out_neighbors(v,g)
		for w in children
			rindex[w] == 0 && visit(w)
			rindex[w] < rindex[v] && begin rindex[v] = rindex[w]; root = false end
		end
		if root
			index -= 1
			while !isempty(S) && rindex[v] <= rindex[S[end]]
				w = pop!(S)
				rindex[w] = c
				index -= 1
			end
			rindex[v] = c
			c -= 1
		else
			push!(S,v)
		end
	end

	for v in vertices(g)
		if rindex[v] == 0
			visit(v)
		end
	end

	return rindex
end

# Pearce algorithm state
type state{T}
	v::T
	stage::T
	root::Bool
end

# Pearce version of Tarjan algorithm - iterative version
function pearce_iterative{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	n = length(vertices(g))
	rindex = zeros(T,n)
	S = T[]
	index = 1
	c = n-1

	function visit(v)
		states = state[]
		current_state = state(v,convert(T,0),true)
		push!(states, current_state)

		@label start
		while !isempty(states)
			current_state = pop!(states)
			active_loop = false

			if current_state.stage == 0
				rindex[current_state.v] = index
				index += 1
				active_loop = true
			end
			children = out_neighbors(current_state.v,g)
			for w in children
				if active_loop
					if rindex[w] == 0
						current_state.stage = w
						push!(states,current_state)

						new_state = state(w,convert(T,0),true)
						push!(states,new_state)
						@goto start
					end
					rindex[w] < rindex[current_state.v] && begin rindex[current_state.v] = rindex[w]; current_state.root = false end
				elseif current_state.stage == w
					active_loop = true
					rindex[w] < rindex[current_state.v] && begin rindex[current_state.v] = rindex[w]; current_state.root = false end
				end
			end
			if current_state.root
				index -= 1
				while !isempty(S) && rindex[current_state.v] <= rindex[S[end]]
					w = pop!(S)
					rindex[w] = c
					index -= 1
				end
				rindex[current_state.v] = c
				c -= 1
			else
				push!(S,current_state.v)
			end
		end
	end

	for v in vertices(g)
		if rindex[v] == 0
			visit(v)
		end
	end

	return rindex
end

# get the main SCC
function get_core{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	sccs = pearce_iterative(g)
	scc_ids = union(sccs,[])
	id_size = Dict{T,T}()

	for id in scc_ids
		id_size[id] = 0
	end

	for id in sccs
		id_size[id] += 1
	end

	sizes = collect(values(id_size))
	sort!(sizes)
	msize = sizes[end]

	# find the max SCC id
	mid = 0
	for id in keys(id_size)
		id_size[id] == msize && begin mid = id; break end
	end

	sids = T[]
	counter = convert(T,1)
	for id in sccs
		id == mid && push!(sids,counter)
		counter = convert(T,counter+1)
	end

	@debug("# core vids: ", length(sids))

	return subgraph(g,sids)
end
# get the main SCC and write it to the specified file (MGS format)
# write the subgraph in MGS format v2
function get_core_streamed{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},sccs::Array{T,1},name::String)
	scc_ids = union(sccs,[])
	id_size = Dict{T,T}()

	for id in scc_ids
		id_size[id] = 0
	end

	for id in sccs
		id_size[id] += 1
	end

	sizes = collect(values(id_size))
	sort!(sizes)
	msize = sizes[end]

	# find the max SCC id
	mid = 0
	for id in keys(id_size)
		id_size[id] == msize && begin mid = id; break end
	end

	sids = T[]
	counter = convert(T,1)
	for id in sccs
		id == mid && push!(sids,counter)
		counter = convert(T,counter+1)
	end

	@debug("# core vids: ", length(sids))

	subgraph_streamed(g,sids,name)
end

# get the reverse graph (same graph with all edge directions reversed)
function get_reverse_graph{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	ng = adjlist(T, is_directed=true)
	vs = vertices(g)

	# same set of vertices
	for v in vs
		add_vertex!(ng,v)
	end

	# inverse the edge directions
	for v in vs
		children = out_neighbors(v,g)
		for c in children
			add_edge!(ng,c,v)
		end
	end

	return ng
end

# get in-degree of g vertices
# returns a dictionary (vertex_id -> in-degree)
function get_vertex_in_degrees{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	vin = Dict{T,T}()
	for v in vertices(g)
		ovs = out_neighbors(v,g)
		for o in ovs
			if !haskey(vin,o)
				vin[o] = convert(T,1)
			else
				vin[o] = convert(T,vin[o]+1)
			end
		end
	end
	return vin
end

# get in- and out- degree arrays of specified graph
function get_in_out_degrees{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	vin = get_vertex_in_degrees(g)
	in_degrees = T[]
	out_degrees = T[]
	for v in vertices(g)
		push!(in_degrees,vin[v])
		push!(out_degrees,convert(T,length(out_neighbors(v,g))))
	end
	return in_degrees,out_degrees
end

# get the avg out-degree of the specified set of visited nodes
#
# p_avg: current average
# nb_steps: number of points used to compute p_avg
#
# NB: to get the avg in-degree of visited nodes, one can use the reverse graph of g
function get_avg_out_degree{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, visited::Array{T,1}, p_avg::Float64=float64(-1), np_steps::Uint64=uint64(0))
	sum = float64(0)
	for v in visited
		sum += length(out_neighbors(v,g))
	end
	if p_avg != -1
		return (p_avg*np_steps+sum)/(np_steps+length(visited))
	else
		return sum/length(visited)
	end
end

# get the list of colinks
function list_colinks{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},filename::String)
	vs = vertices(g)
	n = length(vs)
	io = T[]
	# number of unique colinks per node
	un2 = zeros(Uint64,n)
	count = 0
	total = 0
	ttotal = 0
	bmax_size = 1000
	# delete file
	if isfile(filename)
		rm(filename)
	end
	for v in vs
		# get greater children and parents
		gc = filter(x->x>v,out_neighbors(v,g))
		gp = filter(x->x>v,out_neighbors(v,rg))
		inter = intersect(gc,gp)
		for t in inter
			append!(io,[v,t])
			ttotal += 1
		end
		un2[v] = length(inter)
		total += 1
		count += 1
		if count == bmax_size
			@debug("writing $bmax_size computed vertices [total # colinks: $ttotal | total # of computed vertices: $total (", (total/n*100),"%)]")
			f = open(filename, "a")
			write(f,io)
			io = T[]
			close(f)
			count = 0
		end
	end
	serialize_to_file(un2,"vertex-un2.jld")
	@debug("Total # of unique colinks: ", sum(un2))
	f = open(filename, "a")
	write(f,io)
	close(f)
end

# listing triangles - map
function list_triangles_map{T<:Unsigned}(v::T,g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},dpos::Dict{T,T})
	candidates = (T,T,T)[]
	gc = filter(x->dpos[x]>dpos[v],out_neighbors(v,g))
	gp = filter(x->dpos[x]>dpos[v],out_neighbors(v,rg))
	#@debug("# greater children: ", length(gc))
	#@debug("# greater parents: ", length(gp))
	for c in gc
		# rank(child) > rank(v)
		for p in gp
			# child != parents && rank(parent) > rank(v)
			c != p && push!(candidates,(v,c,p))
		end
	end
	return candidates
end

# listing triangles - reduce
function list_triangles_reduce{T<:Unsigned}(candidates::Array{(T,T,T),1},g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},io::Array{T,1})
	count = convert(T,0)
	for c in candidates
		cc = out_neighbors(convert(T,c[2]),g)
		sort!(cc)
		if binary_search(cc,c[3]) != -1
			append!(io,[c[1],c[2],c[3]])
			count = convert(T,1+count)
		end
	end
	return count
end

# listing triangles - map&reduce
function list_triangles_mapreduce{T<:Unsigned}(v::T,g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},dpos::Dict{T,T},io::Array{T,1})
	count = convert(T,0)
	vpos = dpos[v]
	gc = filter(x->dpos[x]>vpos,out_neighbors(v,g))
	gp = filter(x->dpos[x]>vpos,out_neighbors(v,rg))
	for c in gc
		cc = out_neighbors(c,g)
		sort!(cc)
		for p in gp
			if c != p
				# NB: searchsortedfirst is unreliable!
				#if searchsortedfirst(cc,p) <= length(cc)
				if binary_search(cc,p) != -1
					append!(io,[v,c,p])
					count = convert(T,1+count)
				end
			end
		end
	end
	return count
end

# initialize un2 (# of colink candidates) and un3 (# of triangle candidates) array
function init_un23_arrays{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},dpos::Dict{T,T})
	@debug("initializing # colink candidates + # triangle candidates arrays")
	vs = vertices(g)
	n = length(vs)
	un2 = Uint64[]
	un3 = Uint64[]
	for v in vs
		# total order relation 1
		# get the set of greater children of vertex v
		gc1 = filter(x->x>v,out_neighbors(v,g))
		push!(un2, length(gc1))

		# total order relation 2
		vpos = dpos[v]
		gc2 = filter(x->dpos[x]>vpos,out_neighbors(v,g))
		gp = filter(x->dpos[x]>vpos,out_neighbors(v,rg))
		p = length(gp)
		c = length(gc2)
		i = length(intersect(gp,gc2))
		push!(un3, (p*c-i))
	end
	return un2,un3
end

# initialize mn2 (maximum # of colinks) and mn3 (maximum # of triangles) array
function init_mn23_arrays{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}})
	@debug("initializing max # colinks + max # triangles arrays")
	vs = vertices(g)
	n = length(vs)
	mn2 = Uint64[]
	mn3 = Uint64[]
	for v in vs
		parents = out_neighbors(v,rg)
		children = out_neighbors(v,g)
		p = length(parents)
		c = length(children)
		i = length(intersect(parents,children))
		# maximum number of triangles:
		# = (p_minus_c * c_minus_p) + i * (p_minus_c + c_minus_p + (i - 1))
		# = p*c - i
		push!(mn2, (p+c-i))
		push!(mn3, (p*c-i))
	end
	return mn2,mn3
end

# get the list of triangles
function list_triangles{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},filename::String)
	n = length(vertices(g))
	io = T[]
	# actual number of unique triangles
	un3 = zeros(Uint64,n)
	count = 0
	total = 0
	ttotal = 0
	bmax_size = 10000000
	# delete triangle file if needed
	if isfile(filename)
		rm(filename)
	end
	# maximum number of triangles
	#mn3 = load_serialized("vertex-mn3.jld")
	mn2,mn3 = init_mn23_arrays(g,rg)
	#serialize_to_file(mn2,"vertex-mn2.jld")
	#serialize_to_file(mn3,"vertex-mn3.jld")
	@debug("sorting mn3 array")
	#smn3 = load_serialized("smn3.jld")
	smn3 = bottom_up_sort!(mn3)
	#serialize_to_file(smn3,"smn3.jld")
	@debug("mn3 array sorted")
	@debug("loading position dictionary")
	#dpos = load_serialized("dpos.jld")
	dpos = Dict{T,T}()
	pos = 1
	for v in smn3
		dpos[v] = pos
		pos += 1
	end
	#serialize_to_file(dpos,"dpos.jld")
	@debug("dictionary loaded")
	for v in smn3
		#@debug("Checking vertex: $v")
		#tic()
		#candidates = list_triangles_map(convert(T,v),g,rg,dpos)
		#@debug("# candidates: ", length(candidates))
		#nadded = list_triangles_reduce(candidates,g,io)
		nadded = list_triangles_mapreduce(convert(T,v),g,rg,dpos,io)
		un3[v] = nadded
		#@debug("# found triangles for vertex $v: $nadded")
		#@debug(toc())
		ttotal += nadded
		count += nadded
		total += 1
		if count >= bmax_size
			@debug("writing $count triangles [total # triangles: $ttotal | total # of computed vertices: $total (", (total/n*100),"%)]")
			f = open(filename, "a")
			write(f,io)
			io = T[]
			close(f)
			count = 0
		end
	end
	serialize_to_file(un3,"vertex-un3.jld")
	@debug("Total # of unique triangles: ", sum(un3))
	# write remaining triangles
	f = open(filename, "a")
	write(f,io)
	close(f)
end

# load the colinks distribution from the specified file
function load_colinks_distribution{T<:Unsigned}(size::T,filename::String)
	f = open(filename,"r")
	d = Dict{T,T}()
	while !eof(f)
		v1 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		v2 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		if haskey(d,v1)
			d[v1] += 1
		else
			d[v1] = 1
		end
		if haskey(d,v2)
			d[v2] += 1
		else
			d[v2] = 1
		end
	end
	close(f)
	c_stats = zeros(T,size)
	for i in 1:size
		if haskey(d,i)
			c_stats[i] = d[i]
		else
			c_stats[i] = 0
		end
	end
	return c_stats
end

# load the triangles distribution from the specified file
function load_triangles_distribution{T<:Unsigned}(size::T,filename::String)
	f = open(filename,"r")
	d = Dict{T,T}()
	while !eof(f)
		v1 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		v2 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		v3 = reinterpret(T,read(f,Uint8,sizeof(T)))[1]
		if haskey(d,v1)
			d[v1] += 1
		else
			d[v1] = 1
		end
		if haskey(d,v2)
			d[v2] += 1
		else
			d[v2] = 1
		end
		if haskey(d,v3)
			d[v3] += 1
		else
			d[v3] = 1
		end
	end
	close(f)
	t_stats = zeros(T,size)
	for i in 1:size
		if haskey(d,i)
			t_stats[i] = d[i]
		else
			t_stats[i] = 0
		end
	end
	return t_stats
end

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

# get a ball centered at the specified vertex
function get_forward_ball{T<:Unsigned}(v::T,g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},radius::Int=2,p::Float64=1)
	# vertex ids of the ball
	subids = T[]
	push!(subids,v)
	# set of explored vertices
	explored = T[]
	for t in 1:radius
		so = T[]
		for u in subids
			if !(u in explored)
				append!(so,out_neighbors(u,g))
				push!(explored,u)
			end
		end
		if length(so) > 0
			#  if the probability is 1, add the child
			if p == 1
				append!(subids,unique(so))
			else
				# we are at the first level of the BFS
				if t == 1
					append!(subids,unique(so))
				else
					so = unique(so)
					#added = false
					# select children at random with an exponentially decreasing probability
					for s in so
						if rand() <= p^(t-1)
							push!(subids,s)
							#added = true
						end
					end
					# if no child was added, add at least one at random
					#!added && push!(subids,so[rand(1:length(so))])
				end
			end
		end
	end
	return unique(subids)
end

# get the array of clustering coefficients
function get_clustering_coefficients{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},ntriangles::Array{T,1},density::Float64=-1.)
	vs = vertices(g)
	n = length(vs)
	ccs = zeros(Float64,n)
	for v in vs
		parents = out_neighbors(v,rg)
		children = out_neighbors(v,g)
		p = length(parents)
		c = length(children)
		i = length(intersect(parents,children))
		# maximum number of triangles:
		# = (p_minus_c * c_minus_p) + i * (p_minus_c + c_minus_p + (i - 1))
		# = p*c - i
		npt = p*c-i
		# if there is only 1 colink, no triangle can be formed
		if npt == 0
			ccs[v] = 0
		else
			# correct the result according to the specified density
			if density > 0
				ccs[v] = (ntriangles[v] - density*npt)/(npt*(1-density))
			else
				ccs[v] = ntriangles[v]/npt
			end
		end
	end
	return ccs
end

# get the array of colink coefficients
function get_colink_coefficients{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},rg::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},ncolinks::Array{T,1},density::Float64=-1.)
	vs = vertices(g)
	n = length(vs)
	ccs = zeros(Float64,n)
	for v in vs
		parents = out_neighbors(v,rg)
		children = out_neighbors(v,g)
		p = length(parents)
		c = length(children)
		i = length(intersect(parents,children))
		# maximum number of colinks:
		# = p+c-i
		npc = p+c-i
		# correct the result according to the specified density
		if density > 0
			ccs[v] = (ncolinks[v] - density*npc)/(npc*(1-density))
		else
			ccs[v] = ncolinks[v]/npc
		end
	end
	return ccs
end
