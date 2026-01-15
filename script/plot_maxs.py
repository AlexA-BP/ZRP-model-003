from pathlib import Path
import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import scienceplots
plt.style.use('science')

class Plotting():
    pass

def power_law(t, a, k):
    return a*t**(1/k)

fpath = Path("./data/test_symmetric.h5")

with h5.File(fpath, 'r+') as fid:

    dset = fid["post_processing/maximums"]
    
    if "visualization" not in fid:
        group_id = fid.create_group("visualization") 

    ydata = dset[:]
    xdata = range(len(ydata))
    cutoff = int(1e4)

    popt, pcov = curve_fit(
        power_law, 
        xdata[:cutoff], 
        ydata[:cutoff]
    )
    print(f"a = {popt[0]}, \n k = {popt[1]}")
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.loglog(xdata, ydata)
    ax.loglog(xdata, power_law(xdata, *popt))
    ax.set_title("Maximum size of condensate")
    ax.set_xlabel("time t [a.u.]")
    ax.set_ylabel(r"$\max(\{n_i\})(t)$")
    
    txtstr = (r"fitfunc: $a t^{\frac{1}{k}}$ " + "\n"
              + f"a = {popt[0]:.3f} " + r"$\pm$" 
              + f"{np.sqrt(pcov[0, 0]):.3f}" + "\n"
              + f"k = {popt[1]:.3f} " + r"$\pm$" 
              + f"{np.sqrt(pcov[1, 1]):.3f}" + "\n")
    ax.text(0.05, 0.8, txtstr, 
            transform=ax.transAxes,
            bbox=dict(boxstyle="square",
                      ec="black",
                      fc="white"))

    plt.savefig((fpath.parent)/"max_growth_scaling_symmetric_hopping.png", 
                dpi=300,
                format="png")

    plt.show()



    