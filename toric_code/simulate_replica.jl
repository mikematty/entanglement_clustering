include("lattice.jl")
include("swap_tools.jl")
importall Lattice
using Match

function step!(lat::TCLattice)
  # first have to figure out what kind of move to do
  if((lat.n_spinons > 0) && (rand() > 0.5))
    move_index::Int = rand(1:lat.n_spinons)
    to_move::Index = lat.spinons[move_index]
    moveto::Index = (0,0)
    nns::Vector{Index} = [(1,0),(0,1),(-1,0),(0,-1)]
    shuffle!(nns)
    for dir in nns
      if(move_spinon!(lat,move_index,dir) == 1) return end
    end # for dir in ...
  else 
    flip_plaquette!(lat,rand(1:lat.L),rand(1:lat.L))
  end #if(lat. ...
end # step(...

function swap(lat1::TCLattice,lat2::TCLattice,edges)::Int
  lat1_copy = deepcopy(lat1)
  lat2_copy = deepcopy(lat2)
  for ind in edges
    lat1_copy[ind...],lat2_copy[ind...] = lat2_copy[ind...],lat1_copy[ind...]
  end # for ind in subedges

  if ((count_spinons(lat1_copy) != lat1_copy.n_spinons) ||
      (count_spinons(lat2_copy) != lat2_copy.n_spinons))
    return 0
  else
    return 1
  end
end # function swap(...

function main()
  L::Int = 20
  square_min::Int = 1#div(L,2)
  square_max::Int = 2#L-1
  loop_min::Int = 1
  loop_max::Int = L-2
  rect_min::Int = L-8
  rect_max::Int = L-2
  ncollect::Int = 5000
  n_edge_excitations::Int = 0#40
  nsteps::Int = 1000*ncollect
  nequil::Int = 1000

  lat1 = TCLattice(L)
  for i in 1:L lat1[i,div(L,2),2] = -1 end
  for i in 1:L lat1[div(L,2),i,1] = -1 end
  for i in 1:n_edge_excitations create_random_spinon!(lat1) end
  lat2 = deepcopy(lat1)
  
  for step in 1:nsteps
    step!(lat1)
    step!(lat2)
    if (mod(step,ncollect) == 0)
      rep1, rep2 = make_replicas(lat1),make_replicas(lat2)
      ex1, ex2 = is_spinon_mat(rep1),is_spinon_mat(rep2)
      ul1, ul2 = upper_left_counts(ex1),upper_left_counts(ex2)
      xors = edge_xors(rep1,rep2)
      dirs1, dirs2 = directional_es(ex1,xors),directional_es(ex2,xors)
      res = []
      for lr in 1:L-1, lc in 1:L-1
        avg = 0
        for cr in 1:L, cc in 1:L
          avg += swap_memoized(ul1,ul2,dirs1,dirs2,ex1,ex2,cr,cc,lr,lc)
        end
        append!(res,avg/(L*L))
      end
      println(res)
    end # if (mod(...
  end # for step in ...
end # main(...

main()
