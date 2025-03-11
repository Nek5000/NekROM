import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import re

if len(sys.argv) != 2:
    print("Usage: python plot_coef.py <directory>")
    sys.exit(1)

directory = sys.argv[1]
match = re.search(r'deim_pts(\d+)', directory)
if match:
    deim_pts_number = int(match.group(1))
else:
    print("Error: 'deim_pts' followed by a number not found in the directory string.")
    sys.exit(1)
ucoef_path = os.path.join(directory, "ucoef")


data = np.loadtxt(ucoef_path)
adtr1 = np.reshape(data, (2000, 41), order='F')
t = np.linspace(9640.04,9720, 2000)

data = np.loadtxt("../ops/uk")
uk = np.reshape(data, (2001, 41))
t_snap = np.linspace(9640,9720,2001)

#print(adtr1[:,1])
#print(uk[:,1])

mode = 1
fig1, ax1 = plt.subplots(tight_layout=True)
ax1.plot(t_snap,(uk[:, mode]), 'k-', label='FOM projection')
ax1.plot(t,adtr1[:, mode], 'r--', label='DEIM-ROM')
ax1.legend(loc=0)
ax1.set_xlabel(r"Time $t$")
ax1.set_ylabel(r"ROM coefficient")
ax1.set_title(f"DEIM points:{deim_pts_number}")
fname = f"coef.pdf"
fig1.savefig(os.path.join(directory, fname), format='pdf')
plt.close(fig1)
