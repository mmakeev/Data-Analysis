import re
import numpy as np
from monty.io import zopen
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt

__author__ = 'Rasha Atwi'
__version__ = '0.1'
__email__ = 'rasha.atwi@tufts.edu'
__date__ = 'Aug 8, 2019'


def get_msd(filename):
    """
    This function reads a text file containing msd data and saves it to a dictionary.
    :param filename: name of the file containing time (fs) and msd data (Angstrom^2)
    :return: dictionary with keys of time (s) and values of msd data (cm^2)
    """
    #msd_patt = re.compile(r"(Step)\s*(\w*\s*)(\wmsd)\D")
    #tab_patt = re.compile(r"^\s*(\d+)\s*(\d+.\d*)\s*(\d+.\d*)\s*(\d+.\d*)\D")

    msd = {}

    read_msd = False
    number_of_cols = None
    idx = None
    with zopen(filename) as f:
        for line in f:
            if line.startswith('Step') and 'msd' in line:
                read_msd = True
                line_split = line.split()
                number_of_cols = len(line_split)
                idx = [i for i, j in enumerate(line_split) if 'msd' in j]
                continue
            if read_msd:
                line_split = line.split()
                if len(line_split) != number_of_cols or re.search('(?![eE])[a-zA-Z]', line): #exclude e from string search
                    break
                # Converting time from fs to s and msd from Angstrom^2 to cm^2
                time = int(line_split[0]) / 10 ** 15
                disp = [float(line_split[i]) / 10 ** 16 for i in idx]
                msd[time] = disp
    return msd

def get_diff(filename):
    """
    This function fits the saved msd data with linear regression to calculate diffusion coefficient in cm^2/s.
    :param filename: name of file containing time (fs) and msd (Angstrom^2)
    :return: diffusion coefficient (cm^2/s)
    """
    msd = get_msd(filename)
    time_new = [*msd]
    time_array = np.array(time_new).reshape(-1, 1)
    disp_array = np.transpose(np.array(list(msd.values())))
    for d in disp_array:
        model = LinearRegression().fit(time_array, d)
        slope = model.coef_ /6
        intercept = model.intercept_
        r_sq = model.score(time_array, d)
        print('coefficient of determination:', r_sq)
        print('intercept:', intercept)
        print('diffusion coefficient:', slope)
        plt.plot(np.transpose(time_array)[0], d, color='g')
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), frameon=False, ncol=3)
    plt.xlabel('time (s)')
    plt.ylabel('msd ($cm^2$)')
    plt.savefig('MSD'+'.eps', format='eps', dpi=1000)
    plt.show()
    plt.close()
    #TODO: add name of msd before display

