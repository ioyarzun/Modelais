import matplotlib.pyplot as plt
import pandas as pd

def energy_profile_plot(options, dir, code, data_unopt, data_opt=None):
    """
    This function recieves 1 or 2 datasets and draws the energy profiles.
    Returns a graph of the energy profiles.
    """
    plt.clf()
    # Parameter data_unopt: data of the not optimized structure
    data = pd.read_csv(data_unopt, sep="\s+", header=None)
    line_norm, = plt.plot(data[0], data[1], linewidth=1, label='Not Optimized')

    if options.optimize:
        # Parameter data_opt: data of the optimized structure
        data_o = pd.read_csv(data_opt, sep="\s+", header=None)
        line_opt, = plt.plot(data_o[0], data_o[1], linewidth=1, label='Optimized')
        plt.legend(handles=[line_norm, line_opt])

    plt.title(code + ' energy profile')
    plt.xlabel('Residue')
    plt.ylabel('DOPE_energy')
    plt.axhline(linewidth=1, linestyle=':', color='r')

    # plt.show()
    plt.savefig(dir + "/" + code + '_EnergyProfile_plot.png', dpi=300)

    return
