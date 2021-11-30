import numpy as np, matplotlib.pyplot as plt
import umap, pickle

def load_complex(fname):
  f = open(fname,'r')
  res = []
  for line in filter(lambda x: len(x) > 0, f.readlines()):
    vals = map(lambda x: x.replace("+-","-"), line.split())
    vals = map(lambda x: [complex(x).real,complex(x).imag],vals)
    res += [np.array([val for sublist in vals for val in sublist])]
  #res=[np.array(map(lambda x: complex(x),line.split()),dtype=np.complex128)\
  #     for line in f.readlines()]
  res = filter(lambda x: len(x) > 0, res)
  f.close()
  return np.vstack(res)

def main():
  files = ["../basis_check_outfiles/t2_ex0_data1500_config1.txt",\
           "../basis_check_outfiles/t2_ex10_data1500_config1.txt"]
  legend =["GS", "ES1","ES2","ES3","ES4","ES5"]
  handle = open("basis_check_embedded/t2_ex10.p",'wb')

  swaps = [i for i in map(lambda f: load_complex("output_files/"+f), files)]
  labels = np.hstack([np.array([i]*len(swaps[i])) for i in range(len(swaps))])
  color_labs = [labels == i for i in np.unique(labels)]
  colors = ['r','g','b','y','m','k']
  all_swaps = np.vstack(swaps)
  
  encoded = umap.UMAP(n_neighbors=400).fit_transform(all_swaps)
  pickle.dump((encoded,color_labs),handle)
  handle.close()

  _,ax = plt.subplots()
  for i, color in enumerate(color_labs):
    ax.scatter(encoded[color,0],encoded[color,1],c = colors[i],alpha=0.6,label=legend[i])

  ax.legend(loc=4)
  plt.show()

if __name__=='__main__':main()
