from vmcsim import *
import sys, progressbar, time

def main():
  l,dim = 100, 1 # system size and dimensionality
  ne = l//2 # system filling
  sub_l = 6 # subsystem size for swap calculation
  t1,t2 = 1.0,4.0 # alternating hopping strengths
  nequils,ncollect = 1000,1000
  nsteps = int(1500*ncollect)
  fname = "basis_check_outfiles/t3_ex10_data1500_config1.txt"
  
  # select random excitations (0 for ground state)
  excited = np.random.choice(np.arange(ne),10)
  # Hamiltonian for band insulator
  def onedinsulator(mom,ind):
    k = (4.0*np.pi/l)*(mom+1.0)
    factor = np.exp((ind+2.0-(ind%2))/2.0*k*1.0j)
    nrg = t2+t1*np.exp(k*1j)
    sign = 1.0 if (mom in excited) else -1.0
    return sign*factor if (ind%2) else factor*nrg/np.abs(nrg)
  
  # Hamiltonian for alternate basis in appendix
  def basischeck(mom,ind):
    k = (4.0*np.pi/l)*(mom+1.0)
    factor = np.exp((ind+2.0-(ind%2))/2.0*k*1.0j)
    sign = 1.0 if (mom in excited) else -1.0
    asite = 1j*t2*np.sin(k)
    bsite = t1+t2*np.cos(k)+sign*np.sqrt(t1*t1+t2*t2+2*t1*t2*np.cos(k))
    return factor*bsite if (ind%2) else factor*asite

  sim1 = SimState(l,ne,dim,basischeck)#onedinsulator)
  sim2 = SimState(l,ne,dim,basischeck)#onedinsulator)
  sim2.wvfn = np.copy(sim1.wvfn)
  sim2.rlspc = np.copy(sim1.rlspc)
  #sim2.wvfn_inv = np.linalg.inv(sim2.wvfn)
  #sim2.occ_number = sim1.occ_number[:]
  #sim2.ind_space = sim1.ind_space[:]
  for i in range(nequils): 
    sim1.step()
    sim2.step()
  
  bar = progressbar.ProgressBar(maxval=50, \
      widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
  bar.start()
  
  res = []
  cs = np.zeros(sim1.l,dtype=np.complex128)
  subsystems = [[i%l for i in range(j,j+sub_l)] for j in range(l)]
  ar1,ar2 = 0.0,0.0
  t0 = time.time()
  for nrun in range(nsteps):
    truth1 = sim1.step()
    truth2 = sim2.step()
    ar1 += float(truth1)
    ar2 += float(truth2)
    #####  print progress bar ######
    progress = float(nrun)/nsteps;
    pos = int(round(50*progress))
    bar.update(pos)
    ################################
    if (nrun % ncollect) == 0:
      res += [np.array([swap_op(sim1,sim2,sub) for sub in subsystems],\
              dtype=np.complex128)]
      #for x in xrange(sim1.l): cs[x] += two_point(sim1,x)
  bar.finish()
  tf = time.time()
  print(tf-t0)
  #cs[x] /= nsteps/ncollect
  #for c in cs: print np.abs(c)
  handle = open(fname,'w')#file(fname,'a')
  np.savetxt(fname,np.vstack(res))
  handle.close()
  print(ar1/nsteps*100.0, ar2/nsteps*100.0)

  return 0

if __name__ == '__main__': main()
