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

using DataStructures, Logging

@Logging.configure(level=DEBUG)

##### custom implementation of QuickSort
function swap{T<:Integer}(A::Array{T,1},i::T,j::T)
	t = A[i]
	A[i] = A[j]
	A[j] = t
end

function partition{T<:Integer}(A::Array{T,1},R::Array{T,1},l::T,h::T)
	pvalue = A[h]
	sindex = l
	for j in l:(h-convert(T,1))
		if A[j] < pvalue
			(sindex != j && A[sindex] != A[j]) && begin swap(A,sindex,j); swap(R,sindex,j) end
			sindex += convert(T,1)
		end
	end
	(sindex != h && A[sindex] != A[h]) && begin swap(A,sindex,h); swap(R,sindex,h) end
	return sindex
end

# iterative quicksort
#
# @return permutation array and sort A in ascending order
function quicksort_iterative!{T<:Integer}(A::Array{T,1})
	s = Stack(Tuple{T,T})
	l = convert(T,1)
	h = convert(T,length(A))
	push!(s,(l,h))
	R = collect(l:h)
	while !isempty(s)
		l,h = pop!(s)
		p = partition(A,R,l,h)
		if (p-1) > l
			push!(s,(l,(p-1)))
		end
		if (p+1) < h
			push!(s,((p+1),h))
		end
	end
	return R
end
#####

##### custom implementation of MergeSort

# merge sort (ascending order)
#
# @returns permutation array R
#
# NB: to get the sorted array sA: sA[i] = A[R[i]] 
#     R[i] is thus the index in the original array of element at new index i in the sorted array
function bottom_up_sort{T<:Integer}(A::Array{T,1})
	n = convert(T,length(A))
	B = zeros(T,n)
	# permutation array
	R = collect(convert(T,1):n)
	S = zeros(T,n)
	width = 1
	while width < n
		i = 0
		while i < n
			bottom_up_merge(A,convert(T,i),convert(T,min(i+width,n)),convert(T,min(i+2*width,n)),B,R,S)
			i += 2*width
		end
		A = copy(B)
		R = copy(S)
		width = 2*width
	end
	return R
end

function bottom_up_merge{T<:Integer}(A::Array{T,1},iLeft::T,iRight::T,iEnd::T,B::Array{T,1},R::Array{T,1},S::Array{T,1})
	i0 = iLeft
	i1 = iRight
	for j in (iLeft+1):iEnd
		if i0 < iRight && (i1 >= iEnd || A[i0+1] <= A[i1+1])
			B[j] = A[i0+1]
			S[j] = R[i0+1]
			i0 += 1
		else
			B[j] = A[i1+1]
			S[j] = R[i1+1]
			i1 += 1
		end
	end
end
#####

# custom binary search
#
# search x position in array A
#
# NB: array A is assumed to be sorted in ascending order
function binary_search{T<:Integer}(A::Array{T,1}, x::T)
	low = 1 
	high = length(A)
	while true
		if high == low
			if A[low] == x
				return low
			else
				return -1
			end
		else
			p = low + floor(Int, (high-low)/2)
			if x == A[p]
				return p
			elseif x > A[p]
				if p < high
					low = p+1
				else
					return -1
				end
			else
				if p > low
					high = p-1
				else
					return -1
				end
			end
		end
	end
end

##### custom implementation of Huffman encoding 

# get sorted array
#
# A: initial array
# R: permutation array
function get_sorted_array{T}(A::Array{T,1}, R::Array{T,1}, asc::Bool=true)
	S = zeros(Int,length(A))
	n = length(A)
	# ascending order
	if asc
		for i in 1:n
			S[i] = A[R[i]]
		end
	# decreasing order
	else
		for i in 1:n
			S[n+1-i] = A[R[i]]
		end
	end
	return S
end

# binary tree node type
abstract AbstractNode
type EmptyNodeType <: AbstractNode end
const EmptyNode = EmptyNodeType()

type Node{T} <: AbstractNode
    key::T
    left::AbstractNode
    right::AbstractNode
end

# encode binary tree
#
# the tree is encoded in
# S: bits array
# D: array of leaf node values 
function encode_tree!{T}(root::AbstractNode, S::BitArray{1} , D::Array{T,1})
	if root == EmptyNode
        	push!(S, 0)
        	return
	else
		push!(S, 1)
		push!(D, root.key)
		encode_tree!(root.left,S,D) 
		encode_tree!(root.right,S,D) 
	end
end

# get Huffman prefix codes dictionary
#
# C: dictionary (bitarray -> value::T)
function get_huffman_codes!{T}(root::AbstractNode, C::Dict{BitArray{1},T}, B::BitArray{1})
	if root.key != 0
        	C[B] = root.key
	else
		B1 = copy(B)
		get_huffman_codes!(root.left, C, push!(B1,false)) 
		B2 = copy(B)
		get_huffman_codes!(root.right, C, push!(B2,true)) 
	end
end

# decode binary tree
#
# S: bits array
# D: array of leaf node values 
function decode_tree!{T}(S::BitArray{1}, D::Array{T,1})
	length(S) == 0 && return EmptyNode
	b = shift!(S)
	if b == 1 
        	key = shift!(D)
        	root = Node{T}(key,EmptyNode,EmptyNode)
        	root.left = decode_tree!(S,D)
        	root.right = decode_tree!(S,D)
        	return root
	end
	return EmptyNode
end

# decode values
#
# C: code -> value dictionary
function decode_values{T}(tree::Node{T}, CDATA::BitArray{1})
	children = T[]
	cnode = tree
	for i in 1:length(CDATA)
		if cnode.left == EmptyNode && cnode.right == EmptyNode
			push!(children,cnode.key)
	 		cnode = tree
		end
		if CDATA[i] == 0
			cnode = cnode.left
		elseif CDATA[i] == 1
			cnode = cnode.right
		end
	end
	return children
end

# @return huffman tree
#
# A is assumed to have a length >= 2
# A[i] is the value associated to element having index i (e.g. A[i] could be the in-degree of vertex i)
#
# conventions
# -> lowest child is assigned to left leaf, and highest child to right leaf
# -> 0: left branch, 1: right branch
function huffman_encoding{T<:Unsigned}(A::Array{T,1})
	# get sorted array in increasing order
	# NB: instead of poping elements, we use shift
	## merge sort 
	#R = bottom_up_sort(A)	
	#S = get_sorted_array(A,R)
	## quick sort (A is sorted)
	S = copy(A)
	R = quicksort_iterative!(S)	
	# second queue
	N = Array{Node,1}()
	V = Array{T,1}()
	# creating initial tree
	# NB: lowest node on the left
	lKey = shift!(S)
	rKey = shift!(S)
	# get corresponding elements
	lKey2 = shift!(R)
	rKey2 = shift!(R)
	# create new leaf nodes
	lNode = Node{T}(lKey2,EmptyNode,EmptyNode)
	rNode = Node{T}(rKey2,EmptyNode,EmptyNode)
	nKey = lKey+rKey
	# creating internal node
	nNode = Node{T}(convert(T,0),lNode,rNode)
	# NB: lenght(N) == 0 at the very end of the computation
	while length(S) > 0 || length(N) > 0
		# if some elements remain in the queues, insert latest node to the node queue
		pos = searchsortedfirst(V,nKey)
		insert!(V, pos, nKey)
		insert!(N, pos, nNode)
		# compare lowest values of both queues for node 1
		if length(S) > 0 && S[1] <= V[1]
			key1 = shift!(S)
			key12 = shift!(R)
			# create new leaf node
			node1 = Node{T}(key12,EmptyNode,EmptyNode)
		else
			key1 = shift!(V)
			node1 = shift!(N)
		end
		# compare lowest values of both queues for node 2
		if length(S) > 0 && S[1] <= V[1]
			key2 = shift!(S)
			key22 = shift!(R)
			# create new leaf node
			node2 = Node{T}(key22,EmptyNode,EmptyNode)
		else
			key2 = shift!(V)
			node2 = shift!(N)
		end
		nKey = key1+key2
		if key1 > key2
			nNode = Node{T}(convert(T,0),node2,node1)
		else
			nNode = Node{T}(convert(T,0),node1,node2)
		end
	end
	return nNode
end
#####
