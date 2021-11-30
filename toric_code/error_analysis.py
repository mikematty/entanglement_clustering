import sys, os, pickle, numpy as np, matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Times New Roman'
matplotlib.rcParams['mathtext.it'] = 'Times New Roman:italic'
matplotlib.rcParams['mathtext.bf'] = 'Times New Roman:bold'

CGND = "red"#(0,109/255.,219/255.)
CEXC = "green"#(146/255.,0,0)

def main():
  # User part
  colors = [CGND,CEXC,'b','y','m']
  root = ""#"pickled_embeddings/embd_"
  files=["25l_%dpercent.p"%(i) for i in [5,10,15,20]]#,10,15,20]]
  rows = [0,0,0,0]#,0,1,1]
  cols = [0,1,2,3]#,1,0,1]
  titles = [""]#["Spinon Density %d %%"%(i) for i in [5,10,15,20]]#,10,15,20]]
  #titles = ["One Excitation","Two Excitations","Three Excitations"
  files = ["final_embedded/25l_20nex.p"]

  # Hopefully shouldn't have to touch this part much (except maybe plotting)
  data = []
  for i in range(len(files)):
    handle = open(root+files[i],'rb')
    u = pickle._Unpickler(handle)
    u.encoding = 'latin1'
    data += [u.load()]
    handle.close()

  # The error analysis here only works for two colors
  fig,axs = plt.subplots()#len(np.unique(rows)),len(np.unique(cols)))
  for j,(encoded,color,counts) in enumerate(data):
    ax = axs#[cols[j]]#axs[rows[j],cols[j]]
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title(titles[j],fontname="Times New Roman")
    kmeans = KMeans(n_clusters = 2).fit(encoded)
    ######
    if j in [0,1,2]:
      print(len(color))
      print(sum(1-kmeans.labels_[color[0]]))
      print(sum(kmeans.labels_[color[1]]))
      print(sum(kmeans.labels_))
    ######
    error = 0.0
    for i,color in enumerate(color):
      ax.scatter(encoded[color,0],encoded[color,1],c = colors[i],alpha=0.6)
      cluster = int(round(np.mean(kmeans.labels_[color])))
      cur_counts = counts[color]
      assignments = (kmeans.predict(encoded[color]) != cluster)*cur_counts
      error += float(sum(assignments))
    error = 100.0*(1.0-error/sum(counts))
    ax.scatter(kmeans.cluster_centers_[:,0],kmeans.cluster_centers_[:,1],\
              c='black',marker='D',s=8)
    #if j != 0:
    #ax.set_xlabel("Classification Accuracy %0.2f %%"%(error),fontname="Times New Roman")
    #else: 
    #  ax.set_xlabel("Classification Accuracy 58.88 %",fontname="Times New Roman")
    print(error)

  #plt.suptitle("$L = 20$",fontname="Times New Roman")
  plt.tight_layout()
  plt.savefig("tc_cluster.pdf")
  plt.show()

main()
