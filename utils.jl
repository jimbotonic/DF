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

using DataStructures

##### custom implementation of QuickSort
function swap{T<:Integer}(A::Array{T,1},i::T,j::T)
	t = A[i]
	A[i] = A[j]
	A[j] = t
end

function partition{T<:Integer}(A::Array{T,1},R::Array{T,1},l::T,h::T)
	pvalue = A[h]
	sindex = l
	for j in l:(h-1)
		if A[j] < pvalue
			(sindex != j && A[sindex] != A[j]) && begin swap(A,sindex,j); swap(R,sindex,j) end
			sindex += 1
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
# S: bits array
# D: array of leaf node values 
function encode_tree!{T}(root::Node{T}, S::BitArray{1} , D::Array{T,1})
	if root == EmptyNode
        	push!(S,0)
        	return
	else
		push!(S,1)
		push!(D,root.key)
		encode_tree!(root.left,S,D) 
		encode_tree!(root.right,S,D) 
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

# @return huffman tree
#
# NB: A is assumed to have a length >= 2
function huffman_encoding{T<:Unsigned}(A::Array{T,1})
	# get sorted array in increasing order
	# NB: instead of poping elements, we use unshift
	## merge sort 
	#R = bottom_up_sort(A)	
	#S = get_sorted_array(A,R)
	## quick sort (A is sorted)
	R = quicksort_iterative!(A)	
	S = A
	# second queue
	N = Array{Node}()
	V = Array{T}()
	# initial tree
	# NB: lowest node on the left
	lKey = unshift!(S)
	rKey = unshift!(S)
	lNode = Node{T}(lKey,EmptyNode,EmptyNode)
	rNode = Node{T}(rKey,EmptyNode,EmptyNode)
	nKey = lKey+rKey
	nNode = Node{T}(nKey,lNode,rNode)
	push!(N,nNode)
	push!(V,nKey)
	# NB: lenght(N) == 0 at the very end of the computation
	while length(S) > 0 || length(N) > 0
		# compare lowest values of both queues for node 1
		if length(S) > 0 && S[1] <= V[1]
			key1 = unshift!(S)
			node1 = Node{T}(key1,EmptyNode,EmptyNode)
		else
			key1 = unshift!(V)
			node1 = unshift!(N)
		end
		# compare lowest values of both queues for node 2
		if length(S) > 0 && S[1] <= V[1]
			key2 = unshift!(S)
			node2 = Node{T}(key2,EmptyNode,EmptyNode)
		else
			key2 = unshift!(V)
			node2 = unshift!(N)
		end
		nKey = key1+key2
		if key1 > key2
			nNode = Node{T}(nKey,node2,node1)
		else
			nNode = Node{T}(nKey,node1,node2)
		end
		# if the node list is not empty, insert new node to the node queue
		if length(N) > 0
			pos = searchsortedfirst(V,nKey)
			insert!(V, pos, nKey)
			insert!(N, pos, nNode)
		end
	end
	return nNode
end
#####
