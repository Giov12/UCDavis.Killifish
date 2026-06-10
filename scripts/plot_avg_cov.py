#!/bin/env python3

import matplotlib.pyplot as plt

colors = ["#ff0000", "#ff8700", "#ffd300", "#deff0a", "#a1ff0a", 
         "#0aff99", "#0aefff", "#147df5", "#580aff", "#be0aff"] # vibrant fusion color palette

orgs = ["lpar", "lgoo", "fxen", "fnota", "foli", "fgra", "frat", "fsim", "fcat", "fnot"]
covs = [18.9869, 20.0542, 41.4827, 47.9647, 50.3268, 56.0892, 59.338, 72.1609, 75.005, 78.372]


def make_bar() -> int:
    """make a bar graph"""

    global orgs, covs

    plt.bar(orgs, covs, color="#2E8B57")
    plt.title("Average Coverage")
    plt.ylabel("Coverage (x)")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    plt.savefig("avg_coverage.pdf")

    return 0

def main() -> int:
    """execute the singular function"""

    make_bar()

    return 0

if __name__ == "__main__":
    main()
