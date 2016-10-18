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

# MGS adjacency graph format - version 1
#
# unsigned integers stored in big endian format (byte and bit levels)
# positions: 0-based
#
# graph.index: vid1 (1 word), pos1 (1 word) | ...
# graph.data: vid1, vid2, vid3, ...
#
# NB: vertices do not need to be contiguous
#

# MGS adjacency graph format - version 2
#
# unsigned integers stored in big endian format (byte and bit levels)
# positions: 1-based
#
# graph.index: pos1 (1 word) | pos2 | ...
# graph.data: vid1 (1 word), vid2, vid3, ...
#
# NB: vertex ids need to be contiguous. In case, a vertex has no child, pos[k] = pos[k+1]

using Graphs, DataStructures, HDF5, JLD

# load serialized JLS data
function load_jls_serialized(filename::String)
	x = open(filename, "r") do file
		deserialize(file)
	end
	return x
end

# serialize data to JLS format
function serialize_to_jls(x, filename::String)
	open(filename, "w") do file
		serialize(file, x)
	end
end

# load serialized JLD data
function load_jld_serialized(name::String, filename::String)
  x = jldopen(filename, "r") do file
    read(file, name)
  end
	return x
end

# serialize data to JLS format
function serialize_to_jld(x, name::String, filename::String)
	jldopen(filename, "w") do file
		write(file, name, x)
	end
end

# get the ordered dictionary vid -> startpos
function load_mgs1_graph_index{T<:Unsigned}(ipos::OrderedDict{T,T},filename::String)
	f = open(filename, "r")
	while !eof(f1)
		id = read(f,Uint8,sizeof(T))
		pos = read(f,Uint8,sizeof(T))
		# Julia is 1-based -> +1
		ipos[reinterpret(T,reverse(id))[1]] = reinterpret(T,reverse(pos))[1] + 1
	end
	close(f)
end

# get the set of positions
function load_mgs2_graph_index{T<:Unsigned}(pos::Array{T,1},filename::String)
	f = open(filename, "r")
	while !eof(f)
		p = read(f,Uint8,sizeof(T))
		push!(pos, reinterpret(T,reverse(p))[1])
	end
	close(f)
end

# write index file
function write_mgs1_graph_index{T<:Unsigned}(ipos::OrderedDict{T,T}, filename::String)
	f = open(filename, "w")
	for p in ipos
    # reinterpret pair (vid,pos) in an array of bytes
		bytes = reinterpret(Uint8, [p])
		write(f, reverse(bytes))
	end
	close(f)
end

# write index file
function write_mgs2_graph_index{T<:Unsigned}(pos::Array{T,1}, filename::String)
	f = open(filename, "w")
	for p in pos
    # reinterpret pos in an array of bytes
		bytes = reinterpret(Uint8, [p])
		write(f, reverse(bytes))
	end
	close(f)
end

# get the array of graph children (MGSv1 & MGSv2)
function load_graph_data{T<:Unsigned}(children::Array{T,1},filename::String)
	f = open(filename, "r")
	while !eof(f)
		child = read(f,Uint8,sizeof(T))
		push!(children,reinterpret(T,reverse(child))[1])
	end
	close(f)
end

# write the array of graph children (MGSv1 & MGSv2)
function write_graph_data{T<:Unsigned}(children::Array{T,1}, filename::String)
	f = open(filename, "w")
	for c in children
		bytes = reinterpret(Uint8, [c])
		write(f, reverse(bytes))
	end
	close(f)
end

# load graph in format MGS v1
function load_mgs1_graph{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},name::String)
	ipos = OrderedDict{T,T}()
	children = T[]

	load_mgs1_graph_index(ipos,"$name.index")
	load_graph_data(children,"$name.data")

	# ipos is an ordered dictionary
	ks = collect(keys(ipos))
	# vertex set
	vs = sort(union(children,ks))

	# add vertices
	for i in 1:length(vs)
        	add_vertex!(g,convert(T,i))
	end

	# instantiate dictionary
	oni = Dict{T,T}()

	counter = convert(T,1)
	for v in vs
		oni[v] = counter
		counter += convert(T,1)
	end

	# add edges
	for i in 1:length(ks)
		source = oni[ks[i]]
		# if we reached the last parent vertex
		if i == length(ks)
			pos1 = ipos[ks[i]]
			pos2 = length(children)
		else
			pos1 = ipos[ks[i]]
			pos2 = ipos[ks[i+1]]-1
		end
		for p in pos1:pos2
			target = oni[children[p]]
			add_edge!(g,source,target)
		end
	end
	return g,oni
end

# load graph in format MGS v3
function write_mgs3_graph{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},name::String)
  # 8 bytes string + 4 bytes position: 'MGSv3   ' +  <32bits pos>
  version = 0x4d47537633202020

  vs = vertices(g)

  for v in vs

  end
end

# load graph in format MGS v2
function load_mgs2_graph{T<:Unsigned}(g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}},name::String)
  pos = T[]
  children = T[]

  load_mgs2_graph_index(pos,"$name.index")
	load_graph_data(children,"$name.data")

	#@debug("#pos:",length(pos))
	#@debug("#children:",length(children))

	# vertex set
	vs = range(1,length(pos))

	# add vertices
	for i in 1:length(vs)
        	add_vertex!(g,convert(T,i))
	end

	# add edges
	for i in 1:length(vs)
		source = convert(T,i)
		# if we reached the last parent vertex
		if i == length(vs)
			pos1 = pos[i]
			pos2 = length(children)
		else
			pos1 = pos[i]
			pos2 = pos[i+1]-1
		end
		for p in pos1:pos2
			target = children[p]
			add_edge!(g,source,target)
		end
	end
end

# load triangles from text-formatted file
function load_triangles{T<:Unsigned}(::Type{T},filename::String)
	f = open(filename,"r")
	a = (T,T,T)[]
	while !eof(f)
		te = split(readline(f)[2:(end-2)],',')
		v1 = parseint(T,te[1])
		v2 = parseint(T,te[2])
		v3 = parseint(T,te[3])
		push!(a,(v1,v2,v3))
	end
	close(f)
	return a
end
