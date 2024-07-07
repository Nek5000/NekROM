import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.unicode']=True
import matplotlib.pyplot as plt

# number of basis function
nb = 60
# number of basis function in fom
nbs = 60
# number of snapshots
K = 500

data = np.zeros((nb,K))
data_fom = np.zeros((nbs,K))
for i in range(500):
   filename=str(i+1).zfill(4)+".out"

# load the data with NumPy function loadtxt
   data[:,i] = np.loadtxt(filename)
   data_fom[:,i] = np.loadtxt("../500snaps_200nb_fom/"+filename)

mean = np.zeros((nb,))
var = np.zeros((nb,))
mean_fom = np.zeros((nb,))
var_fom = np.zeros((nb,))

for i in range(nb):
   mean[i] = np.sum(data[i,:])/len(data[0])
   var[i] = np.sum((data[i,:]-mean[i])**2)/(len(data[0])-1)
for i in range(nb):
   mean_fom[i] = np.sum(data_fom[i,:])/len(data_fom[0])
   var_fom[i] = np.sum((data_fom[i,:]-mean_fom[i])**2)/(len(data_fom[0])-1)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(1)
plt.plot(1+np.arange(nb),mean,'ro',markersize=5)
plt.plot(1+np.arange(nb),mean_fom,'b-o',lw=1,markersize=5,markerfacecolor="None")
plt.xlabel(r'\textbf{n}')
plt.ylabel(r'$\langle a_n \rangle_s$',fontsize=16)
plt.title(r"Sample mean of coefficient $\{a^{j}_n\}_j$",fontsize=16) 
axes = plt.gca()
ymin = -5
ymax = 20
axes.set_ylim([ymin,ymax])
plt.savefig('mean.png')

plt.figure(2)
plt.plot(1+np.arange(nb),var,'ro',markersize=5)
plt.plot(1+np.arange(nb),var_fom,'b-o',lw=1,markersize=5,markerfacecolor="None")
plt.xlabel(r'\textbf{n}')
plt.ylabel(r'$\mathcal{V}_s(a_n)$',fontsize=16)
plt.title(r"Sample variance of coefficient $\{a^{j}_n\}_j$",fontsize=16) 
axes = plt.gca()
ymax = 80
axes.set_ylim([ymin,None])
plt.savefig('variance.png')
   
plt.show()

