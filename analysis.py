import os, sys, argparse, math
import pandas as pd

R = 1.987204e-3 # Gas Constant in kcal/mol
T = 310
PKI2DG = -1.3637

def to_dG(pKi):
    """
    Given a pKi, compute a delta G

    Parameters:
        pKi, the calculated binding affinity, -log10(Ki)
    Returns:
        dG at 310 K
    """
    return -1.3637 * float(pKi)

# Read and analyze the results of AP-Net-dG and/or MMGBSA
if __name__ == '__main__':

    ## ==> Read in the arguments <== ##
    parser = argparse.ArgumentParser(description='A module to analyze the results of different binding energy prediction softwares')
    parser.add_argument('setname', help='[string] The name of the validation set you are evaluating the different methods on (ex: 4mxo_site_1_struct_1_S985_val')
    parser.add_argument('--skip_mmgbsa', help='[bool] Skip MMGBSA in the analysis?', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])
    setname = args.setname
    do_mmgbsa = not args.skip_mmgbsa

    dg_path = os.environ.get('APNETDG')
    if dg_path is None: raise Exception("Environment variable $APNETDG is not set.")

    results = {
        'system' : [],
        'ref_dG' : [],
        'apnet_dG' : []
    }

    if do_mmgbsa:
        results['mmgbsa_dG'] = []

    dg_data = pd.read_pickle(f'{dg_path}/datasets/{setname}/dimers.pkl')
    for n in range(len(dg_data)):
        results['system'].append(dg_data['system'][n])
        results['ref_dG'].append(to_dG(dg_data['label'][n]))

    dg_pred = pd.read_csv(f'apnetdg/{setname}/preds.csv')
    for n in range(len(dg_pred)):
        results['apnet_dG'].append(to_dG(dg_pred['AP_net_dG_pKi'][n]))

    if do_mmgbsa:
        mm_res = pd.read_csv(f'mmgbsa/{setname}/pocket-out.csv')
        for n in range(len(mm_res)):
            results['mmgbsa_dG'].append(mm_res['r_psp_MMGBSA_dG_Bind'][n])

    df = pd.DataFrame(data=results, dtype='object')
    df_dir = os.path.join('analysis', setname)
    if not os.path.isdir(df_dir): os.makedirs(df_dir)
    df.to_csv(f'{df_dir}/results.csv', index=False)