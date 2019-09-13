import sys, os
sys.path.append(os.environ['RICHESTRMCL'])
import numpy as np
import configparser as cp
import generate_mass_catalog as gmc
import draw_true_properties as dtp
import draw_observed_properties as dop
import hmf_library as hmflib
import choose_model as cm
import mock_setup as ms


def print_highlight(text):
    print("#############################")
    print(text)
    print("#############################")
    return


def main():
    # Load the config file in the above directory
    # print_highlight("MANUALLY SETTING CONFIG FILE")
    cfg_name = sys.argv[1]#'../config/' + sys.argv[1]
    print_highlight(cfg_name)
    cfgin = cp.ConfigParser()
    cfgin.read(cfg_name)

    # Set up the number of runs
    total_runs = int(sys.argv[2])
    start_at = int(sys.argv[3])

    # Load lookup tables and other constants
    sinfo, cosmology = ms.import_cosmology(cfgin)
    splines = ms.build_splines(cfgin, sinfo)
    mgas_rescale = 1.e12

    # Set up object to draw lamobs
    drawer = dop.Converter()
    drawer.setup_iCDF_grid()

    for run_number in np.arange(total_runs):
        run_number += start_at
        print("Run %i/%i" % (run_number+1, total_runs+start_at))

        print("> Selecting model")
        file_names = ms.build_file_names(cfgin, run_number)
        model = cm.draw_model(cfgin, sinfo, file_names, splines)

        print("> Drawing masses")
        gmc.generate_mass_catalog(cfgin, model, sinfo, file_names)

        print("> Drawing true properties")
        model['rich_pivot'] = np.log(cfgin['General'].getfloat('rich_pivot'))
        model['mgas_pivot'] = np.log(cfgin['General'].getfloat('mgas_pivot'))
        catalog = dtp.draw_true_properties(cfgin, model, sinfo, file_names, splines)

        print("> Drawing observed properties")
        model['stdobs_lnmgas'] = 0.14
        catalog = dop.draw_observed_properties(catalog, model, drawer, mgas_rescale)
        ms.save_observed_catalog(catalog, file_names)

        print("> Finding the richest N")
        dtp.find_richest_n(catalog, cfgin, file_names)


    return





if __name__ == "__main__":
    main()
