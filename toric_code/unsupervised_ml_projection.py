import numpy as np, matplotlib.pyplot as plt
import umap, pickle

def main():
  # Files to load for clustering
  files = ["../final_outfiles/t2_ex0_data1000_config1.txt",\
           "../final_outfiles/t2_ex2_data1500_config1.txt",\
           "../final_outfiles/t2_ex2_data1500_config2.txt",\
           "../final_outfiles/t2_ex2_data1500_config3.txt",\
           "../final_outfiles/t2_ex2_data1500_config4.txt",\
           "../final_outfiles/t2_ex2_data1500_config5.txt"]
  legend =["GS1", "ES1","ES2","ES3","ES4","ES5"]
  handle = open("final_embedded/t2_ex2_multiconfig.p",'wb')

  swaps = [i for i in map(lambda f: load_complex("output_files/"+f), files)]
  labels = np.hstack([np.array([i]*len(swaps[i])) for i in range(len(swaps))])
  color_labs = [labels == i for i in np.unique(labels)]
  colors = ['r','g','b','y','m','k']
  all_swaps = np.vstack(swaps)
  
  # Do UMAP projection
  encoded = umap.UMAP(n_neighbors=400).fit_transform(all_swaps)
  pickle.dump((encoded,color_labs),handle)
  handle.close()

  _,ax = plt.subplots()
  for i, color in enumerate(color_labs):
    ax.scatter(encoded[color,0],encoded[color,1],c = colors[i],alpha=0.6,label=legend[i])

  ax.legend(loc=4)
  plt.show()

if __name__=='__main__':main()
