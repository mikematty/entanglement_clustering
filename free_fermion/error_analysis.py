import sys, os, pickle, numpy as np, matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Times New Roman'
matplotlib.rcParams['mathtext.it'] = 'Times New Roman:italic'
matplotlib.rcParams['mathtext.bf'] = 'Times New Roman:bold'

def main():
  # User part
  colors = ['r','g','b','y','m','k']
  root = "basis_check_embedded/"#"output_files/encoded_"
  ex_files = ["t2_ex4.p","t2_ex6.p","t2_ex8.p","t2_ex10.p"]
  #gap_files= ["t0_ex10.p","t1_ex10.p","t2_ex10.p","t3_ex10.p","t4_ex10.p"]
  files = ex_files#+gap_files
  rows = [0,0,0,0,0,1,1,1,1,1]
  cols = [0,1,2,3,4,0,1,2,3,4]
  titles = [""]*10#["$N_{ex}/L = %d$"%(i)+"%" for i in [1,2,3]]+\
           #["$\Delta E = {}$".format(i) for i in [0.0, 0.1, 0.5]]
  #titles = ["One Excitation","Two Excitations","Three Excitations"
  
  # Hopefully shouldn't have to touch this part much (except maybe plotting)
  data = []
  for i in range(len(files)):
    handle = open(root+files[i],'rb')
    u = pickle._Unpickler(handle)
    u.encoding = 'latin1'
    data += [u.load()]
    ############# Added for MS edits ###############
    #data[i] = (data[i][0][int(((data[i][0]).shape[0])*.7):,:],data[i][1])
    ##################################################################
    handle.close()

  # The error analysis here only works for two colors
  fig,axs = plt.subplots()#len(np.unique(rows)),len(np.unique(cols)))
  for j,(encoded,color) in enumerate([data[2]]):#enumerate(data):
    ax = axs#axs[rows[j],cols[j]]
    ax.set_yticks([])
    ax.set_xticks([])
    #ax.set_title(titles[j],fontname="Times New Roman")
    kmeans = KMeans(n_clusters = 2).fit(encoded)
    error = 0.0
    for i,color in enumerate(color):
      ax.scatter(encoded[color,0],encoded[color,1],c = colors[i],alpha=0.6,s=6)
      cluster = int(round(np.mean(kmeans.labels_[color])))
      error += float(sum(kmeans.predict(encoded[color]) != cluster))
    error = 100.0*(1.0-error/len(encoded))
    #ax.scatter(kmeans.cluster_centers_[:,0],kmeans.cluster_centers_[:,1],\
    #          c='black',marker='D')
    print(error)
    #ax.set_xlabel("Classification Accuracy %0.2f %%"%(error),fontname="Times New Roman")

    # Special Plotting Stuff I May have to mess with
    #if(rows[j] == 0) and (cols[j] == 0):
    #  ax.set_ylabel("$\Delta E = 1.0$",fontsize=15)
    #if(rows[j] == 1) and (cols[j] == 0):
    #  ax.set_ylabel("$N_{ex}/L = 10$"+"%",fontsize=15,fontname="Times New Roman")
  #plt.suptitle("$L = 100$, 1D, Two Band Tight Binding Model",fontname="Times New Roman")
  plt.savefig("band_insulator_clusters.pdf")
  plt.show()

main()
