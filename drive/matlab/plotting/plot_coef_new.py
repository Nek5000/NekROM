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

my_case = "cyl"
#my_case = "ldc"

if my_case == "cyl":
    nsteps = int(10*1.25000E+05);#20000; 
    dt     = 4.000000E-03;#0.001;
    iostep = 5000;#250;%500;%10;

    #nsteps = int(1e5)*10;#1.25000E+05;%20000;
    sample_cutoff_time = 4e-3*1.25e05; 
    #dt     = 1.000000E-03;#0.001;
    #iostep = 1000*5#40;#250;%500;%10;
    nb=21
    title="Flow past a cylinder"

    directories = ["../cyl_nb20_results_GROM",
               "../cyl_nb20_results_ndeim_pts256__GROM",
               '../cyl_nb20_results_ndeim_pts256__TR',
               '../cyl_nb20_results_ndeim_pts256_clsdeim_GROM',
               #"nb20_results_deim_pts2000_EFR",
               #"nb20_results_ndeim_pts200__GROM",
               #"nb20_results_ndeim_pts200__CROM", 
               #"nb20_results_ndeim_pts200_clsdeim_GROM",
               #"nb20_results_ndeim_pts200_mclsdeim_GROM",
               #"nb20_results_deim_pts2000_EFR"
               #"nb40_results_sopt_pts200_EFR",
               #"nb40_results_deim_pts2000_GROM", 
               #"nb40_results_deim_pts2000_EFR", 
               #"nb40_results_deim_pts2000_TR", 
               #"nb40_results_deim_pts2000_Leray", 
               #"nb40_results_GROM", "nb40_results_EFR", "nb40_results_TR", "nb40_results_Leray"
               ]


elif my_case == "ldc":
    nsteps = int(10*1e5);#80000;%1.25000E+05;%20000; 
    dt     = 1.000000E-03;#0.001;
    iostep = 5*1000;#500;%250;%500;%10;
    nb = 21
    title = "Lid-driven cavity"

    directories = ["../ldc_nb20_results_GROM",
               "../ldc_nb20_results_ndeim_pts256__GROM",
               '../ldc_nb20_results_ndeim_pts256__TR',
               '../ldc_nb20_results_ndeim_pts256_clsdeim_GROM',
               #"nb20_results_deim_pts2000_EFR",
               #"nb20_results_ndeim_pts200__GROM",
               #"nb20_results_ndeim_pts200__CROM", 
               #"nb20_results_ndeim_pts200_clsdeim_GROM",
               #"nb20_results_ndeim_pts200_mclsdeim_GROM",
               #"nb20_results_deim_pts2000_EFR"
               #"nb40_results_sopt_pts200_EFR",
               #"nb40_results_deim_pts2000_GROM", 
               #"nb40_results_deim_pts2000_EFR", 
               #"nb40_results_deim_pts2000_TR", 
               #"nb40_results_deim_pts2000_Leray", 
               #"nb40_results_GROM", "nb40_results_EFR", "nb40_results_TR", "nb40_results_Leray"
               ]



#nsnapshots = 80000//iostep;
n_io_steps = nsteps//iostep;
io_dt = iostep*dt;

deim_pts_number = 256

labels = ["GROM",
          "GROM-DEIM",
          #"DEIM-EFR", 
          "GROM-DEIM-TR", 
          #"DEIM-Leray", 
          #"DEIM-GROM", 
          #"DEIM-CROM", 
          #"EFR",
          #"TR", 
          "GROM-CLSDEIM",
          #"MCLDEIM-GROM",
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
    ax1.set_xlabel(r"Time $t$")
    ax1.set_ylabel(r"Value of first ROM coefficient")
    #ax1.set_title(f"DEIM points:{deim_pts_number}")
    ax1.set_title(title)

ax1.vlines(100, -3, 3,  label="Sample window cutoff", color='k')
ax1.legend(loc=0)
#data = np.loadtxt("../ops/uk")
#uk = np.reshape(data, (2001, nb))
#ax1.plot(np.arange(0,2001)*io_dt,(uk[:, mode]), 'k-', label='FOM projection')

fname = my_case+ "_" + f"coef.pdf"
#plt.show()
fig1.savefig(os.path.join('./', fname), format='pdf')
#plt.close(fig1)
