"""
This file contains the utilities for plotting the simulation outputs and
reading data from solver output files.
"""
import matplotlib.pyplot as plt


def simple_plot(times, values, selected_species, xlog=False):
    """
    times ... list of timestamps
    values ... dictionary of lists of the same dimension as times
    selected_species ... keys from values to be plotted
    xlog ... specifies if logarithmic x-axis

    You should call plt.show() or plt.savefig(...) after running this method. Optionally, this call can be preceded
    by other plot settings: such as plt.title(...).
    """
    plt.yscale("log")  # set the y-axis in the plot to logarithmic scale
    if xlog: plt.xscale("log")
    plt.xlabel("time [s]")
    plt.ylabel("number of particles")
    for species in selected_species:
        plt.plot(times, values[species], label=species)
    plt.legend()


def plot_with_EN(times, values, selected_species, ENs=None, xlog=False, ylim=None):
    """
    times ... list of timestamps
    values ... dictionary of lists of the same dimension as times
    selected_species ... keys from values to be plotted
    Ens ... list of values of E/N at timestamps from times (must have the same dimension)
    xlog ... specifies if logarithmic x-axis
    ylim ... if specified, sets ylim for the values y-axis

    You should call plt.show() or plt.savefig(...) after running this method. Optionally, this call can be preceded
    by other plot settings: such as plt.title(...).
    """
    fig, ax1 = plt.subplots()

    if xlog: ax1.set_xscale("log")
    ax2 = ax1.twinx()

    ax1.set_yscale('log')
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel("species concentration [m$^{-3}$]")
    ax2.set_ylabel('EN', color='b')

    for species in selected_species:
        ax1.plot(times, values[species], label=species)

    if not ENs: ENs = values['EN']
    ax2.plot(times, ENs, 'b', label="EN")
    ax1.legend()
    ax2.legend()
    if ylim: ax1.set_ylim(ylim)


def read_outfile(filename):
    """
    Reads an output file created by a method from solver.py and translates it into times, values in the same format
    as the solvers from solver.py would return.
    """
    f = open(filename, 'r')
    values = {}
    for line in f.readlines():
        spl = line.split(', ')
        for item in spl:
            item = item.strip()
            name, val = item.split(': ')
            if name not in values:
                values[name] = []
            values[name].append(float(val))
    f.close()

    # time has been parsed as a parameter of values
    times = values['time']
    del values['time']

    return times, values
