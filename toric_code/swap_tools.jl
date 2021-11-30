# matrix in site space where each entry is the number of spinons up and to
# the left of that index 
function upper_left_counts(spinon_mat)
  L = size(spinon_mat)[1]
  ul = zeros(Int,(L,L))
  for row in 1:L
    row_sum = 0
    for col in 1:L
      row_sum += spinon_mat[row,col]
      elt = row_sum
      if (row > 1) elt += ul[row-1,col] end
      ul[row,col] = elt
    end
  end
  return ul
end

function is_spinon_mat(rep)
  L = size(rep)[1]
  epsilon::Array{Int,2} = [is_spinon(rep,i,j) for i in 1:L, j in 1:L]
  return epsilon
end

# array same shape as lattices where each entry is 0 if the two edges are the
# same, and 1 if they are different
function edge_xors(rep1,rep2)
  L = size(rep1)[1]
  x::Array{Int,3} = [1-(rep1[i,j,k] == rep2[i,j,k]) for i in 1:L, j in 1:L, k in 1:2]
  return x
end

# These are used to find the change in the number of spinons on the boundary
# of the subsystem. Takes matrix of whether a site is a spinon and 
# the edge x-or matrix
function directional_es(sp,x)
  e_u,e_d = zeros(Int,size(sp)),zeros(Int,size(sp))
  e_l,e_r = zeros(Int,size(sp)),zeros(Int,size(sp))
  l = size(sp)[1]
  for i in 1:l, j in 1:l
    e_u[i,j] = sum(x[((i==1)?l:(i-1)),jp,2] * (1-2*sp[i,jp]) for jp in 1:j)
    e_d[i,j] = sum(x[i,jp,2] * (1-2*sp[i,jp]) for jp in 1:j)
    e_l[i,j] = sum(x[ip,((j==1)?l:(j-1)),1] * (1-2*sp[ip,j]) for ip in 1:i)
    e_r[i,j] = sum(x[ip,j,1] * (1-2*sp[ip,j]) for ip in 1:i)
  end
  return [e_u,e_d,e_l,e_r]
end

# takes matrix of upper left counts, corner row index, corner column index,
# length in rows of rectangle, length in cols of rectangle
function rectangle_count(ul,cr,cc,lr,lc)
  out::Int = ul[cr+lr,cc+lc]
  if(cr == 1)
    if(cc == 1) return out end
    out -= ul[cr+lr,cc-1]
  elseif(cc == 1)
    out -= ul[cr-1,cc+lc]
  else out = out - ul[cr+lr,cc-1] - ul[cr-1,cc+lc] + ul[cr-1,cc-1]
  end
  return out
end

# takes the directional e things ordered as shown
# takes cr,cc,lr,lc for SUBSYSTEM, NOT expanded boudary
function boundary_count(e_u,e_d,e_l,e_r,cr,cc,lr,lc)
  l = div(size(e_u)[1],2)
  E3Bd = E3Bu = E3Bl = E3Br = 0
  if(cr == 1) cr += l end
  if(cc == 1) cc += l end
  if(lr < l-1)
    E3Bu = e_d[cr-1,cc+lc]-e_d[cr-1,cc-1]
    E3Bd = e_u[cr+lr+1,cc+lc] - e_u[cr+lr+1,cc-1]
  end
  if(lc < l-1)
    E3Bl = e_r[cr+lr,cc-1] - e_r[cr-1,cc-1]
    E3Br = e_l[cr+lr,cc+lc+1] - e_l[cr-1,cc+lc+1]
  end
  res = E3Bd + E3Bu + E3Bl + E3Br
  return res
end

# takes upper left counts, directional boundary deltas, and 
# spinon matrices, and subystem corner row, col and row length and col length
function swap_memoized(ul1,ul2,dir_es1,dir_es2,smat1,smat2,cr,cc,lr,lc)
  r1 = rectangle_count(ul1,cr,cc,lr,lc)
  r2 = rectangle_count(ul2,cr,cc,lr,lc)
  dN2 = r1 - r2 + boundary_count(dir_es2...,cr,cc,lr,lc)
  dN1 = r2 - r1 + boundary_count(dir_es1...,cr,cc,lr,lc)
  return Int((dN1 == 0)*(dN2 == 0))
end

function make_replicas(lat)
  L = lat.L
  replica = TCLattice(2*L)
  replica[1:L,1:L,:] = lat.lat
  replica[(L+1):2*L,(L+1):2*L,:] = lat.lat
  replica[1:L,(L+1):2*L,:] = lat.lat
  replica[(L+1):2*L,1:L,:] = lat.lat
  return replica
end
