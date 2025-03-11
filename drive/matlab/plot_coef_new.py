import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import re
"""
if len(sys.argv) != 2:
    print("Usage: python plot_coef_new.py <directory>")
    sys.exit(1)

directory = sys.argv[1]
match = re.search(r'deim_pts(\d+)', directory)
if match:
    deim_pts_number = int(match.group(1))
else:
    print("Error: 'deim_pts' followed by a number not found in the directory string.")
    sys.exit(1)
"""
nsteps = 100000;#1.25000E+05;%20000; 
dt     = 1.000000E-03;#0.001;
iostep = 500#40;#250;%500;%10;
nb=21

#nsnapshots = 80000//iostep;
n_io_steps = nsteps//iostep;
io_dt = iostep*dt;

deim_pts_number = 200


directories = ["nb20_results_GROM",
               "nb20_results_ndeim_pts200_clsdeim_GROM",
               #"nb20_results_deim_pts2000_EFR"
               #"nb40_results_sopt_pts200_EFR",
               #"nb40_results_deim_pts2000_GROM", 
               #"nb40_results_deim_pts2000_EFR", 
               #"nb40_results_deim_pts2000_TR", 
               #"nb40_results_deim_pts2000_Leray", 
               #"nb40_results_deim_pts2000_CROM", 
               #"nb40_results_GROM", "nb40_results_EFR", "nb40_results_TR", "nb40_results_Leray"
               ]
labels = ["GROM",
          "CLSDEIM-GROM",
          #"DEIM-EFR", 
          #"DEIM-TR", 
          #"DEIM-Leray", 
          #"DEIM-CROM", 
          #"GROM", 
          #"EFR",
          #"TR", 
          #"Leray"
         ]

fig1, ax1 = plt.subplots(tight_layout=True)

t = np.arange(0,n_io_steps)*io_dt

#ax1.set_ylim([-1.5,1.5])
for label, directory in zip(labels, directories):
    ucoef_path = os.path.join(directory, "ucoef")

    print(ucoef_path)

    data = np.loadtxt(ucoef_path)
    adtr1 = np.reshape(data, (n_io_steps, nb), order='F')

    #print(adtr1)
    #exit()
    #t = np.linspace(9640.04,9720, 2000)


    #t_snap = np.linspace(9640,9720,2001)

    #print(adtr1[:,1])
    #print(uk[:,1])

    mode = 1
    ax1.plot(t, adtr1[:, mode], '--', label=label)
    ax1.legend(loc=0)
    ax1.set_xlabel(r"Time $t$")
    ax1.set_ylabel(r"ROM coefficient")
    ax1.set_title(f"DEIM points:{deim_pts_number}")

#data = np.loadtxt("../ops/uk")
#uk = np.reshape(data, (2001, nb))
#ax1.plot(np.arange(0,2001)*io_dt,(uk[:, mode]), 'k-', label='FOM projection')

fname = f"coef.pdf"
plt.show()
#fig1.savefig(os.path.join(directory, fname), format='pdf')
#plt.close(fig1)
