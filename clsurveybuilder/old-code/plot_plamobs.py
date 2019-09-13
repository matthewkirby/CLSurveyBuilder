import matplotlib
import configparser as cp
#matplotlib.use('Agg')
import numpy as np
from pandas import DataFrame
import draw_observed_properties as dop
import mock_setup as ms
from astropy.table import Table
import matplotlib.pyplot as plt
import time
import sys


def make_histograms(ltarr, ztruearr, n_obj, invcdf_grid, maxr):
    conv = dop.Converter()
    for lamtrue, ztrue in zip(ltarr, ztruearr):
        t0 = time.time()
        cat = test_obs_richness(lamtrue, ztrue, invcdf_grid, maxr, n_obj)
        t1 = time.time()
        plt.hist(cat['obs_richness'], bins=200, log=True, normed=True, histtype='step', 
                 label=r'$\lambda_{{true}}={}, z={}$'.format(lamtrue, ztrue))

        # Plot the analytic curve
        lamobs_analytic1 = np.linspace(-10, 3*lamtrue, 10001)
        plt.semilogy(lamobs_analytic1, conv.p_lambda_obs_true(lamtrue, lamobs_analytic1, ztrue), '-k')
        # print("Time for {} is {}".format(lamtrue, t1-t0))

    plt.xlabel(r'$\lambda_{obs}$')
    plt.legend()
    ax = plt.gca()
    ax.set_xlim([0, 300])
    ax.set_ylim([1e-6, 6e-2])
    plt.show()
    return


# def plot_cdfs(ltarr, invcdf_grid):
#     index0, Delta, Deltar = 8.0, 0.1, 0.0001
#     cdf = np.linspace(0.0, 1.0, 10001)
#     for i in range(len(ltarr)):
#         cdf_index = np.floor((ltarr[i] - index0)/Delta).astype(int)
#         lamobs_list = invcdf_grid[cdf_index]
#         plt.plot(lamobs_list, cdf, label=r'$\lambda_{{true}}={}$'.format(ltarr[i]))
#
#     plt.legend()
#     plt.xlabel(r'$\lambda_{obs}$')
#     plt.ylabel(r'CDF')
#     plt.show()
#     return


def test_obs_richness(ltrue, ztrue, invcdf_grid, maxr, n_objs=10000000):
    cfgin = cp.ConfigParser()
    cfgin.read('../mcmc.ini')

    data = np.ones((n_objs, 2))
    data[:,0] = ltrue
    data[:,1] = ztrue
    catalog = DataFrame(data, columns=['true_richness', 'z'])
    catalog = dop.draw_observed_richness(catalog, invcdf_grid, maxr)
    return catalog


def main():
    n_obj = 100000#00
    #ltarr = [23, 55, 74.32, 95]
    ltarr = [25, 50, 75, 95]
    print("Plotting for lamtrue={}".format(ltarr))
    ztrue = [0.1, 0.15, 0.2, 0.25]
    make_histograms(ltarr, ztrue, n_obj, None, None)#invcdf_grid, maxr)
    #plot_cdfs(ltarr, invcdf_grid)



if __name__ == "__main__":
    main()
