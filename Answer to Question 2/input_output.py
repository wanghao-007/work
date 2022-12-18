from __future__ import print_function
from builtins import input
from builtins import str
from builtins import range

import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as AA
import numpy as np


def initial_setup(config_file=None, interactive=True):
    # INPUT DATA AND CONVERT TO NUMPY ARRAYS 'time' and 'abs'
    print()
    print("Data will be written to a file having the same name as the data file, ")
    print("but with _Analysis.txt at the end.")
    name = input("Enter the filename of the data file (it should end in .txt): ")
    print()
    print("If the number of fits selected is 1, then the program will just fit")
    print("the original data set and won't do an extrapolation.")
    steps = input("Enter the number of extrapolation fits (1-10): ")
    steps = int(steps)
    print()
    print("If the enzyme concentration is not known, entering 1 will mean that")
    print("the reported kcat/Km is actually k(obs).")
    enzconc = input("Enter the enzyme concentration in M (e.g 3.2E-9): ")
    enzconc = float(enzconc)
    return name, steps, enzconc


def make_plots(dt):
    # SET UP PLOTS OF DATA
    # NOT EXACTLY SURE ABOUT THE NEXT 4 LINES, But putting data on plot doesn't work without ax defined.
    f = plt.figure()
    f.subplots_adjust(right=0.85)
    ax = AA.Subplot(f, 1, 1, 1)
    f.add_subplot(ax)

    plot_title = dt.name[:-4] + str(dt.loopcount + 1) + "  Relative [S] = " + str(dt.substrate[dt.loopcount])
    x_axislabel = "Time (sec)"
    y_axislabel = "Absorbance"
    plt.title(plot_title)
    plt.xlabel(x_axislabel)
    plt.ylabel(y_axislabel)
    plt.axis([dt.time[0], dt.time[dt.current_num_data - 1], min(dt.abs), min(dt.abs) + 1.1 * (max(dt.abs) - min(dt.abs))])
    plt.plot(dt.time, dt.abs, color='orange')  # The actual data for this fit

    #  The calculated fits superimposed
    pcalc = np.zeros(dt.current_num_data)
    pcalc = dt.init_abs + (dt.vi) * (1 - np.exp(-1. * dt.k * dt.time)) / dt.k  # Method 1 fit
    plt.plot(dt.time, pcalc, 'k', linewidth=1.0)
    pcalc = dt.init_abs2 + (dt.deltaAf) * (1 - np.exp(-1. * dt.k2 * dt.time))  # Method 2 fit
    plt.plot(dt.time, pcalc, 'b', linewidth=1.0)

    # PRINT DATA ON THE PLOT from Fit Method 1
    show_data = "Initial Abs = {:.3e} +/- {:.2e} AU".format(dt.init_abs, np.sqrt(dt.pcov1[0, 0]))
    plt.text(0.6, 0.6, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Initial Rate = {:.3e}  +/- {:.2e} AU/sec".format(dt.vi, np.sqrt(dt.pcov1[1, 1]))
    plt.text(0.6, 0.56, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Rate constant, k = {:.3e} +/- {:.2e} /sec".format(dt.k, np.sqrt(dt.pcov1[2, 2]))
    plt.text(0.6, 0.52, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "[Enzyme] = {:.2e} M".format(dt.enzconc)
    plt.text(0.6, 0.48, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "kcat/Km = {:.3e} +/- {:.2e} /M/sec".format(dt.kcat_Km[dt.loopcount], dt.kcat_Km[dt.loopcount] * np.sqrt(dt.pcov1[2, 2]) / dt.k)
    plt.text(0.6, 0.44, show_data, ha='center', va='center', transform=ax.transAxes)

    # PRINT DATA ON THE PLOT -- From Fit Method 2
    show_data = "Initial Abs = {:.3e} +/- {:.2e} AU".format(dt.init_abs2, np.sqrt(dt.pcov[0, 0]))
    plt.text(0.6, 0.36, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "DeltaAbs = {:.3e}  +/- {:.2e} AU".format(dt.deltaAf, np.sqrt(dt.pcov[1, 1]))
    plt.text(0.6, 0.32, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Final Abs = {:.3e} AU".format(dt.init_abs2 + dt.deltaAf)
    plt.text(0.6, 0.28, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Rate constant, k = {:.3e} +/- {:.2e} /sec".format(dt.k2, np.sqrt(dt.pcov[2, 2]))
    plt.text(0.6, 0.24, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "kcat/Km = {:.3e} +/- {:.2e} /M/sec".format(dt.kcat_Km2[dt.loopcount], dt.kcat_Km2[dt.loopcount] * np.sqrt(dt.pcov[2, 2]) / dt.k)
    plt.text(0.6, 0.20, show_data, ha='center', va='center', transform=ax.transAxes)

    # SAVE THE PLOT AS A PDF FILE
    plfilenm = dt.name[:-4] + "plot" + str(dt.loopcount + 1) + ".pdf"
    plt.savefig(plfilenm, format='pdf')
    print("Saved as ", plfilenm)
    plt.show()
    plt.close(plfilenm)
    return plfilenm


def final_plot(dt):
    fig = plt.figure()
    ax = AA.Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax)

    if dt.kval == 0:
        plot_title = dt.name[:-4] + " -- Extrapolation plot, based on last 5% of data = final absorbance"
    else:
        plot_title = dt.name[:-4] + " -- Extrapolation plot, based on fit final absorbance"
    x_axislabel = "Relative [Substrate]"
    y_axislabel = "Km/kcat (s*mM)"

    plt.plot(dt.x, dt.y, 'ko')  # Putting the data points on the plot
    plt.title(plot_title)
    plt.xlabel(x_axislabel)
    plt.ylabel(y_axislabel)
    plt.axis([0, 1.1 * max(dt.x), 0.9 * min(dt.y), 1.1 * max(dt.y)])
    ty = []
    tx = []
    tx.append(0.)
    ty.append(dt.a)
    tx.append(1.1 * max(dt.x))
    ty.append(dt.a + 1.1 * max(dt.x) * dt.b)
    plt.plot(tx, ty, 'k', linewidth=1.5)  # Putting the line on the plot

    show_data = "Slope = {:0.3e} +/- {:0.1e} /s/mM".format(dt.b, dt.sigb)
    plt.text(0.5, 0.92, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "y-intercept = Km/kcat = {:0.3e} +/- {:0.1e} s mM".format(dt.a, dt.siga)
    plt.text(0.5, 0.88, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "1/intercept = kcat/Km = {:0.3e} +/- {:0.1e} /s/mM".format(1.0 / dt.a, 1.0 / dt.a * dt.siga / dt.a)
    plt.text(0.5, 0.84, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Slope/intercept = {:0.3e}".format(dt.b / dt.a)
    plt.text(0.5, 0.80, show_data, ha='center', va='center', transform=ax.transAxes)

    plfilenm = dt.name[:-4] + "_ExtrapPlot" + str(dt.kval) + ".pdf"
    plt.savefig(plfilenm, format='pdf')
    plt.show()
    plt.close(plfilenm)
    return plfilenm