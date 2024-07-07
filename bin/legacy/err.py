import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm

matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.unicode']=True


colors = cm.rainbow(np.linspace(0, 1, 11))
color_ctr = 0

# number of snapshots
K =100 

# number of basis function you want to show
nb = np.arange(10,60,10)


# load snap_coef
tmp2 = np.loadtxt('../100snap_analysis_stoke_h10/coef_snap')
coef_snap = np.reshape(tmp2,(101,100),order='F')


# set current group parameters
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

count =0
err_H1_mean = np.zeros((len(nb),))
for j in nb:
    # load rom_coef
    tmp1 = np.loadtxt('../100snap_rom_analysis_'+str(j)+'nb_stoke_h10/coef')
    coef_rom = np.reshape(tmp1,(j,100),order='F')
    
    err_H1 = np.zeros((K,))
    for i in range(K):
       err_H1[i] = np.sum((coef_snap[1:j+1,i]-coef_rom[:,i])**2)

    err_H1_mean[count] = np.sqrt(np.sum(err_H1)/K) 
    count = count + 1
    err_H1 = np.sqrt(err_H1)
    
    plt.figure(j)
    plt.semilogy(1+np.arange(K),err_H1,'ro',markersize=5)
    plt.title(r"Error between snap and rom in H1 norm",fontsize=16) 
    plt.legend(['POD-Gal','FOM'],fontsize=16)
    
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.xlabel(r'\textbf{k}',fontsize=20)
    
    ax = plt.figure(j).gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xticks(fontsize=14); plt.yticks(fontsize=14)
    axes = plt.gca()
    ymin = 1e-2; ymax = 1e2 
    xmin = 0; xmax = K+1
    axes.set_ylim([ymin,ymax]); axes.set_xlim([xmin,xmax])

    plt.savefig('err_H1_'+str(j).zfill(2)+'nb.png')
    plt.close()

    if (j == 20 or j == 30 or j == 40 or j == 50 or j == 100):
        plt.figure(0)
        plt.semilogy(1+np.arange(K),err_H1,label=str(j)+r'$nb$',\
                     color=colors[color_ctr], linestyle='-', linewidth=2.5, marker='o')
        color_ctr += 1

plt.legend()
plt.title(r'Error $||\mathcal{P}\hat{u}^k - \tilde{u}^k||_{H^1_0}$ versus $nb$ ',fontsize=14)
plt.xlabel(' kth snapshot ',fontsize=14)
axes = plt.gca()
ymin = 1e-0; ymax = 1e2 
axes.set_ylim([ymin,ymax])
plt.savefig('err_H1.png')

plt.figure()
plt.semilogy(nb,err_H1_mean,'b-o')
axes = plt.gca()
ymin = 1e-0; ymax = 1e2 
axes.set_ylim([ymin,ymax])
plt.title(r'Mean error $\frac{1}{K}\sum^K_{k=1} ||\mathcal{P}\hat{u}^k - \tilde{u}^k\||_{H^1_0}$ ',fontsize=14)
plt.xlabel('-- nb --',fontsize=14)
plt.savefig('mean_err.png')


