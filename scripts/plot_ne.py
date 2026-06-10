#!/bin/env python3

import argparse
import os
import sys
import matplotlib.pyplot as plt
from   glob import glob
from   math import log

input_dir  = ''
psmc_files = list()

class Data:
    def __init__(self) -> None:
        self.time      = list()
        self.pop_sizes = list()


def get_arguments() -> int:
    """get the single argument"""

    d = "Create a line plot of the bootstrapped NE results from PSMC"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--dir", help="directory containing .txt intermediate files from PSMC", required=True)
    args   = parser.parse_args()

    assert os.path.isdir(args.dir), f"Could not locate {args.dir}"

    global input_dir

    input_dir = args.dir

    return 0

def load_files() -> int:
    """collect the text files in the provided path"""

    global input_dir, psmc_files

    if (input_dir[-1] == '/'):
        pattern = f"{input_dir}*.txt"
    else:
        pattern = f"{input_dir}/*.txt"

    psmc_files = glob(pattern)

    if (len(psmc_files) == 0):
        msg = f"No *.txt files were found in {input_dir}"
        sys.exit(msg)

    return 0

def parse_text(text_file: str, data: Data) -> int:
    """parse the text file for the time and effective population size"""

    fh  = open(text_file, 'r')
    msg = f"malformed file detected in {os.path.basename(text_file)}"

    for line in fh:
        if (len(line) == 0 or line[0] == '#'): 
            continue
        fields = line.split('\t')
        assert len(fields) == 5, msg
        time    = float(fields[0])
        time    = 1.0 if time == 0.0 else time # avoid log of 0 when plotting
        popsize = float(fields[1])
        data.time.append(time)
        data.pop_sizes.append(popsize)

    fh.close()

    return 0

def is_bootstrap(text_file: str) -> bool:
    """determine if this is the main psmc file or a bootstrap"""

    bname  = os.path.basename(text_file)
    fields = bname.split('.')
    return fields[-2] != '0'

def plot_ne() -> int:
    """plot the estimates"""

    global psmc_files

    x_label = "Years (g = 1, " + "$\mu$" + " = 0.1 x " + "$\mathrm{10}^{-8}$" + ')'
    y_label = "$N_e$" + " (" + "$\mathrm{10}^{4}$" + ')'

    fig, axe = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
    axe.set_title('')
    axe.set_xlabel(x_label, fontsize = 14)
    axe.set_ylabel(y_label, fontsize = 14)
    axe.set_xlim(10000.0, 100000000.0)
    axe.set_xscale("log")
    axe.set_ylim(0, 750)

    for psmc_file in psmc_files:
        data = Data()
        parse_text(psmc_file, data)
        if (is_bootstrap(psmc_file)):
            alpha = 0.20
            lnwd  = 0.20
        else:
            alpha = 1.0
            lnwd  = 0.75
        axe.step(data.time, data.pop_sizes, color="#FF2800", alpha=alpha, linestyle = '-', linewidth = lnwd)

    plt.savefig("Ne.plot.pdf")

    return 0

def main() -> int:
    """entry point to this little subprogram"""

    # get the input path
    get_arguments()

    # load the input files
    load_files()

    # plot the results
    plot_ne()

    return 0

if __name__ == "__main__":
    main()
