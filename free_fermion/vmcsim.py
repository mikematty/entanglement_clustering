import numpy as np

class SimState(object):
  def __init__(self,l,ne,d,inds_to_elt):
    self.l, self.ne, self.d = l, ne, d
    self.wvfn = np.array([[inds_to_elt(i,j) for j in range(int(l))]\
                         for i in range(ne)],dtype=np.complex128,order='F')
    self.rlspc = np.random.permutation(l)
    self.wvfn = self.wvfn[:,self.rlspc]
    self.wvfn_inv = np.array(np.linalg.inv(self.wvfn[:,:self.ne]),order='C')

  def step(self):
    e_to_move = np.random.randint(0,self.ne)
    moveto = np.random.randint(self.ne,self.l)
    new_cols = [i if i != e_to_move else moveto for i in range(self.ne)]
    r = det_ratio(self.wvfn_inv,self.wvfn[:,new_cols],e_to_move)
    #r=np.linalg.det(self.wvfn[:,new_cols])/np.linalg.det(self.wvfn[:,:self.ne])
    if min(1.0, np.power(np.abs(r),2.0)) > np.random.random():
      cols = [i for i in range(self.l)]
      cols[e_to_move],cols[moveto] = cols[moveto],cols[e_to_move]
      self.wvfn = self.wvfn[:,cols]
      self.rlspc = self.rlspc[cols]
      update_inverse(self.wvfn_inv,self.wvfn,r,e_to_move)
      return True
    return False

def det_ratio(wvfn_i_inv,wvfn_f,col):
  return sum([wvfn_i_inv[col][row]*wvfn_f[row][col] for row in range(wvfn_f.shape[0])])

def update_inverse(wvfn_i_inv, wvfn_f, r, col):
  size = wvfn_f.shape[0]
  for i in range(size):
    if i!= col:
      beta_i = sum([wvfn_i_inv[i,k]*wvfn_f[k,col] for k in range(size)])/r
      wvfn_i_inv[i,:] = wvfn_i_inv[i,:]-beta_i*wvfn_i_inv[col,:]
  wvfn_i_inv[col,:] = wvfn_i_inv[col,:]/r

def two_point(sim,r):
  res = 0.0+0.0j
  for x in range(sim.l):
    y = (x+r)%sim.l
    if ((not x in sim.rlspc[:sim.ne])or(y == x)) and (y in sim.rlspc[:sim.ne]):
      ycol,xcol = np.where(sim.rlspc == y)[0][0],np.where(sim.rlspc == x)[0][0]
      cols = [i if i != ycol else xcol for i in range(sim.ne)]
      res += np.linalg.det(sim.wvfn[:,cols])/np.linalg.det(sim.wvfn[:,:sim.ne])
    return res/sim.l

def swap_op(sim1,sim2,sub):
  ind1 = set(filter(lambda site: site in sim1.rlspc[:sim1.ne], sub))
  ind2 = set(filter(lambda site: site in sim2.rlspc[:sim2.ne], sub))
  ind1,ind2 = ind1-ind2, ind2-ind1
  if not ((len(ind1) == len(ind2)) and (len(ind1) > 0)): return 0.0+0.0j

  res = 1.0+0.0j
  rl1,rl2 = np.copy(sim1.rlspc),np.copy(sim2.rlspc)
  wvfn1,wvfn2 = np.copy(sim1.wvfn),np.copy(sim2.wvfn)
  inv1,inv2 = np.copy(sim1.wvfn_inv),np.copy(sim2.wvfn_inv)
  #old1,old2 = np.copy(wvfn1),np.copy(wvfn2)
  for (st1,st2) in zip(ind1,ind2):
    col11,col12 = np.where(rl1 == st1)[0][0], np.where(rl1 == st2)[0][0]
    col21,col22 = np.where(rl2 == st1)[0][0], np.where(rl2 == st2)[0][0]
    cols1,cols2 = [i for i in range(sim1.l)], [j for j in range(sim2.l)]
    cols1[col11],cols1[col12] = cols1[col12], cols1[col11]
    cols2[col21],cols2[col22] = cols2[col22], cols2[col21]

    rl1, rl2 = rl1[cols1], rl2[cols2]
    wvfn1,wvfn2 = wvfn1[:,cols1], wvfn2[:,cols2]
    r1,r2 = det_ratio(inv1,wvfn1[:,:sim1.ne],col11),det_ratio(inv2,wvfn2[:,:sim2.ne],col22)
    res *= r1*r2

    update_inverse(inv1,wvfn1,r1,col11)
    update_inverse(inv2,wvfn2,r2,col22)

    #res *= np.linalg.det(wvfn1[:,:sim1.ne])*np.linalg.det(wvfn2[:,:sim2.ne])/\
    #      (np.linalg.det(old1[:,:sim1.ne])*np.linalg.det(old2[:,:sim2.ne]))
    #old1,old2 = np.copy(wvfn1),np.copy(wvfn2)
  return res
