import matplotlib.pyplot as plt


def plot(times, values, selected_species):
    plt.yscale("log")  # set the y-axis in the plot to logarithmic scale
    plt.xlabel("time [s]")
    plt.ylabel("number of particles (log scale)")
    for species in selected_species:
        plt.plot(times, values[species], label=species)
    plt.legend()
    plt.show()