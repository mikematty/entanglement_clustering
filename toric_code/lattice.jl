
module Lattice


using Match
export AbstractLattice, SquareLattice, ToricCode, TCLattice, Index, Plaquette, 
       Vedges, get_plaquette, flip_plaquette!, get_vedges, flip_edge!,
       is_spinon, count_spinons, SubLattice, create_random_spinon!, vis,
       move_spinon!, make_rectangular_sublattice, sublattice_union,
       in_bounds, sublattice_edges, count_swap_spinons, edges_and_ownership


# Index type for vertices and edges
Index = Tuple{Int,Int}
# Type for four edges around a plaquette
Plaquette = Tuple{Int,Int,Int,Int}
# Type for four edges incident on a vertex
Vedges = Tuple{Int,Int,Int,Int}
# Type to specify a subsystem
SubLattice = Vector{Index} #Tuple{Index,Index}

# Abstract type for general lattice
abstract AbstractLattice
# Abstract type for 3d square lattice
abstract SquareLattice <: AbstractLattice
abstract ToricCode <: SquareLattice

type TCLattice <: ToricCode
  #Linear Dimension
  L::Int
  # Lattice of spin values on bonds
  # lat[i,j,1] corresponds to the edge to the right of the vertex
  # at (i,j), lat[i,j,2] corresponds to the edge down from vertex (i,j)
  lat::Array{Int,3}
  spinons::Vector{Index}
  n_spinons::Int
  TCLattice(L::Int) = new(L,[1 for i in 1:L, j in 1:L, k  in 1:2],[],0)
end

#type TCSwap <: ToricCode
#  L::Int
#  lat1::Array{Int,3}
#  lat2::Array{Int,3}
#  sub::Nullable{Vector{Int}}
#  TCSwap(L::Int) = new(L,TCLattice(L),TCLattice(L),Nullable{Vector{Int}}())
#end

function Base.size(lat::TCLattice)
  return size(lat.lat)
end

function Base.deepcopy(lat::TCLattice)
  new_lat = TCLattice(lat.L)
  new_lat.lat = deepcopy(lat.lat)
  new_lat.spinons = deepcopy(lat.spinons)
  new_lat.n_spinons = lat.n_spinons
  return new_lat
end

# Vectorize addition for indices
import Base.+
function +(x::Index,y::Index)
  return (x[1]+y[1],x[2]+y[2])
end

# Make indexing happen modularly
function Base.getindex(l::TCLattice,i::Int,j::Int,k::Int)
  return l.lat[mod(i-1,l.L)+1,mod(j-1,l.L)+1,k]
end

function Base.getindex(l::TCLattice,i,j,k)
  if(typeof(i) == Int) i = mod(i-1,l.L)+1 end
  if(typeof(j) == Int) j = mod(j-1,l.L)+1 end
  return l.lat[i,j,k]
end

# Make setting an index happen modularly
function Base.setindex!(l::TCLattice,v::Int,i::Int,j::Int,k::Int)
  l.lat[mod(i-1,l.L)+1,mod(j-1,l.L)+1,k] = v
end

function Base.setindex!(l::TCLattice,v,i,j,k)
  if(typeof(i) == Int) i = mod(i-1,l.L)+1 end
  if(typeof(j) == Int) j = mod(j-1,l.L)+1 end
  l.lat[i,j,k] = v
end

# Bring a modular index in bounds of the lattice size
function in_bounds(ind::Index,l::Int)::Index
  return (mod(ind[1]-1,l)+1,mod(ind[2]-1,l)+1)
end

# index of a plaquette corresponds to the vertex that is at the
# top left corner of the plaquette
# need right and down of top left, down of right, and right of down
# returned as ccw starting from top edge of plaquette
function get_plaquette(lat::ToricCode,i::Int,j::Int)::Plaquette
  return (lat[i,j,1],lat[i,j,2],lat[i+1,j,1],lat[i,j+1,2])
end # get_plaquette(... 

# flip values on bonds surrounding plaquette in place
function flip_plaquette!(lat::TCLattice,i::Int,j::Int)
  for (ii,jj,kk) in [(i,j,1),(i,j,2),(i+1,j,1),(i,j+1,2)]#get_plaquette(lat,i,j)
    lat[ii,jj,kk] *= -1
  end # for (ii,...
end # flip_plaquette!(...

# get values of bonds surrounding vertext at i,j
# need right and down of vert, down of up, and right of left
# returned cw starting from right of vert at i,j
function get_vedges(lat::ToricCode,i::Int,j::Int)::Vedges
  return (lat[i,j,1],lat[i,j,2],lat[i,j-1,1],lat[i-1,j,2])
end # get_vedges(...

# flip edge between sites (i1,j1) and (i2,j2) 
function flip_edge!(lat::TCLattice,i1::Int,j1::Int,i2::Int,j2::Int)
  #println(i2-i1,' ',j2-j1,' ',-(lat.L-1))
  @match (i2-i1,j2-j1) begin
    (0,-1)::Index         => lat[i2,j2,1] *= -1
    (0,1)::Index          => lat[i1,j1,1] *= -1
    (-1,0)::Index         => lat[i2,j2,2] *= -1
    (1,0)::Index          => lat[i1,j1,2] *= -1
    (0,-(lat.L-1))::Index => lat[i1,j1,1] *= -1
    (0,lat.L-1)::Index    => lat[i2,j2,1] *= -1
    (-(lat.L-1),0)::Index => lat[i1,j1,2] *= -1
    (lat.L-1,0)::Index    => lat[i2,j2,2] *= -1 
    _                     => error("ambiguous edge, bad vertices given")
  end # @match...
  # not sure if this is the best place to do this but for now we can
  # update spinon counts when flipping an edge cause that's the only time
  # it should change
  lat.n_spinons = count_spinons(lat)
end # flip_edge!(...

# dispatch for index type arguments
function flip_edge!(lat::TCLattice,i::Index,j::Index)
  flip_edge!(lat,i[1],i[2],j[1],j[2])
end

function move_spinon!(lat::TCLattice,move_index::Int,dir::Index)::Int
  to_move::Index = lat.spinons[move_index]
  if((to_move+dir) in lat.spinons) 
    return 0
  else 
    flip_edge!(lat,to_move,to_move+dir)
    lat.spinons[move_index] = to_move+dir
    return 1
  end # if((to_move...
end # move_spinon!(...

# create a pair of excitations by flipping a random edge
function create_random_spinon!(lat::TCLattice)
  site1::Index = (0,0)
  dir::Index = (0,0)
  sites::Vector{Index} = [(i,j) for i in 1:lat.L for j in 1:lat.L]
  shuffle!(sites)
  nns::Vector{Index} = [(1,0),(0,1),(-1,0),(0,-1)]
  for site in sites
    if (!(site in lat.spinons))
      shuffle!(nns)
      for dir in nns
        if(!(site+dir in lat.spinons))
          # spinon counting is happening in flip_edge! right now. 
          # this probably isn't the best move but it's what's happening now
          flip_edge!(lat,site,site+dir)
          append!(lat.spinons,[site,site+dir])
          return
        end # if(!(site+dir ...
      end # for dir in nns
    end # if (!(site in ...
  end # for site in sites
  error("no viable edge for to create excitations...")
  return
end # create_random_spinon!(...

# determine if a vertex is the endpoint of a string, i.e. a spinon
# returns 1 if it is, 0 if it isn't
function is_spinon(lat::ToricCode,i::Int,j::Int)::Int
  return (prod(get_vedges(lat,i,j)) == 1) ? 0 : 1
  #return (mod(abs(sum(get_vedges(lat,i,j))),4) == 0) ? 0 : 1
end # is_endpoint(...

function is_swap_spinon(lat1::TCLattice,lat2::TCLattice,ownership::Array{Bool,3},
                        i::Int,j::Int)
  L = lat1.L
  nns = [[i,j,1],[i,j,2],[i,mod(j-2,L)+1,1],[mod(i-2,L)+1,j,2]]
  prod1::Int,prod2::Int = 1,1
  for nn in nns
    if ownership[nn...]
      prod1 *= lat1[nn...]
      prod2 *= lat2[nn...]
    else
      prod2 *= lat1[nn...]
      prod1 *= lat2[nn...]
    end
  end
  return ((prod1 == 1) ? 0 : 1), ((prod2 == 1) ? 0 : 1)
end

# count number of spinons in a lattice
function count_spinons(lat::ToricCode)::Int
  spinon_count = sum([is_spinon(lat,i,j) for i in 1:lat.L for j in 1:lat.L])
  if(!iseven(spinon_count)) error("have to have even number of spinons...") end
  return spinon_count
end # count_spinons(...

function count_swap_spinons(lat1::TCLattice,lat2::TCLattice,ownership::Array{Bool,3})
  count1::Int = 0 
  count2::Int = 0
  for i in 1:lat1.L, j in 1:lat1.L
    (is_spinon1,is_spinon2) = is_swap_spinon(lat1,lat2,ownership,i,j)
    count1 += is_spinon1
    count2 += is_spinon2
  end
  return count1, count2
end

function vis(lat::TCLattice)
  L::Int = lat.L
  table::Vector{String} = []
  conv = x -> (x == 1) ? "\033[94m↑\033[0m" : "\033[91m↓\033[0m"
  for row in 1:(2*L+1)
    if(row == 1)
      append!(table,[reduce(*,["   "*conv(lat.lat[L,i,2]) for i in 1:L])*"   "])
    elseif(mod(row,2) == 1)
      append!(table,[reduce(*,["   "*conv(lat.lat[div(row,2),i,2])for i in 1:L])*"   "])
    else
      str::String = ""
      for col in 1:(4*L+3)
        if(mod(col,2) == 1)
          str *= "-"
        elseif(col == 2)
          str *= conv(lat.lat[div(row,2),L,1])
        elseif(mod(col,4) == 2)
          str *= conv(lat.lat[div(row,2),div(col,4),1])
        else
          str *= (is_spinon(lat,div(row,2),div(col,4)) == 1) ? "\033[92m•\033[0m" : "•"
        end # if(mod(col,...
      end # for col in ...
      append!(table,[str])
    end # if(mod(row,...
  end # for row in...
  for line in table
    println(line)
  end
end # function vis(...

# generate the SubLattice, i.e. array of indices
function make_rectangular_sublattice(lx::Int,ly::Int,l::Int,corner::Index)::SubLattice
  res::SubLattice = [corner+(i,j) for i in 0:lx for j in 0:ly]
  for i in 1:length(res), j in 1:2
    res[i] = in_bounds(res[i],l)
  end
  return res
end # function make_rectangular_sublattice

function sublattice_union(s1::SubLattice,s2::SubLattice)::SubLattice
  res::SubLattice = Vector{Index}()
  for s in s1 if(!(s in res)) append!(res,[s]) end end
  for s in s2 if(!(s in res)) append!(res,[s]) end end
  return res
end

# Return sublattice corresponding to s1 \ s2
function sublattice_set_difference(s1::SubLattice,s2::SubLattice)::SubLattice
  res::SubLattice = Vector{Index}()
  for site in s1
    if !(site in s2) append!(res,[site]) end
  end # for site in s1
end

function sublattice_edges(sub::SubLattice,L::Int)
  edges = Vector{}()
  for ind in sub
    #append!(edges,[[ind[1],ind[2],1]])
    #append!(edges,[[ind[1],ind[2],2]])
    append!(edges,[[ind[1],ind[2],:]])
    if(!((in_bounds(ind+(0,-1),L)) in sub))
      append!(edges,[[ind[1],ind[2]-1,1]])
    end
    if(!((in_bounds(ind+(-1,0),L)) in sub))
      append!(edges,[[ind[1]-1,ind[2],2]])
    end
  end # for ind in sub
  return edges
end

function edges_and_ownership(L::Int,subedges)::Array{Bool,3}
  out = Array{Bool,3}(L,L,2)
  for i in 1:L, j in 1:L, k in 1:2
    if [i,j,k] in subedges
      out[i,j,k] = false
    else 
      out[i,j,k] = true
    end
  end
  return out
end

end # Module Lattice
