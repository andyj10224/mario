import os, sys, argparse
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import schrodinger.structure as structure
import schrodinger.structutils.smiles as smiles

R = 1.987204e-3 # Gas Constant in kcal/mol
T = 310
PKI2DG = -1.3637

def to_dG(pKi):
    """
    Given a pKi, compute a delta G

    Parameters:
        pKi [float] : the calculated binding affinity, -log10(Ki)
    Returns:
        dG at 310 K
    """
    return -1.3637 * float(pKi)

def write_summary(name, errs, data, exp, file):
    """
    A helper function to print summary statistics to a file between different binding energy predictions.

    Parameters:
        name (str) : Name of the scoring method (i.e. Docking, AP-Net-dG, MMGBSA)
        errs (list[float]) : The errors for each data point
        data (list[float]) : The predicted dG of the scoring method
        exp (list[float]) : The experimental dG
        file (file, NOT str) : The file object to write the summary to
    """
    me = np.mean(errs)
    rmse = np.sqrt(np.mean(np.square(errs)))
    mae = np.mean(np.abs(errs))

    sr = spearmanr(data, exp)[0]
    r = np.corrcoef(data, exp)[0,1]
    r2 = r**2

    file.write(f'\t{name.upper()} STATISTICS\t => ME: {me:8.3f}, RMSE: {rmse:7.3f}, MAE: {mae:7.3f}, R: {r:6.3f}, R^2: {r2:5.3f}, SPEARMANR: {sr:6.3f}\n')

# Read and analyze the results of AP-Net-dG and/or MMGBSA
if __name__ == '__main__':

    ## ==> Read in the arguments <== ##
    parser = argparse.ArgumentParser(description='A module to analyze the results of different binding energy prediction softwares')
    parser.add_argument('setname', help='[string] The name of the validation set you are evaluating the different methods on (ex: S32_val)')
    parser.add_argument('--skip_mmgbsa', help='[bool] Skip MMGBSA in the analysis?', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])
    setname = args.setname
    do_mmgbsa = not args.skip_mmgbsa

    dg_path = os.environ.get('APNETDG')
    if dg_path is None: raise Exception("Environment variable $APNETDG is not set.")

    results = {
        'system' : [],
        'smiles' : [],
        'exp_dG' : [],
        'exp_std' : [],
        'dock_dG' : [],
        'dock_err' : [],
        'apnet_dG' : [],
        'apnet_err' : []
    }

    if do_mmgbsa:
        results['mmgbsa_dG'] = []
        results['mmgbsa_err'] = []

    dg_data = pd.read_pickle(f'{dg_path}/datasets/{setname}/dimers.pkl')
    ndata = len(dg_data)
    for n in range(ndata):
        results['system'].append(dg_data['system'][n])
        results['exp_dG'].append(to_dG(dg_data['label'][n]))

    prime_datapath = f'docking/{setname}/pocket_pv.maegz'
    for n, st in enumerate(structure.StructureReader(prime_datapath)):
        if n == 0: continue
        dock_dG = st.property['r_i_glide_gscore']
        smiles_gen = smiles.SmilesGenerator()
        smiles_str = smiles_gen.canonicalize(smiles_gen.getSmiles(st))
        results['smiles'].append(smiles_str)
        if 'r_sd_std_pKi' in st.property.keys():
            std_exp = st.property['r_sd_std_pKi']
        else:
            std_exp = st.property['i_sd_std_pKi']
        results['exp_std'].append(abs(PKI2DG) * std_exp)
        results['dock_dG'].append(dock_dG)
        results['dock_err'].append(dock_dG - results['exp_dG'][n-1])

    dg_pred = pd.read_csv(f'apnetdg/{setname}/preds.csv')
    for n in range(ndata):
        pred_dG = to_dG(dg_pred['AP_net_dG_pKi'][n])
        exp_dG = results['exp_dG'][n]
        results['apnet_dG'].append(pred_dG)
        results['apnet_err'].append(pred_dG - exp_dG)

    if do_mmgbsa:
        mm_res = pd.read_csv(f'mmgbsa/{setname}/output.csv')
        for n in range(len(mm_res)):
            mm_dG = mm_res['r_psp_MMGBSA_dG_Bind'][n]
            exp_dG = results['exp_dG'][n]
            results['mmgbsa_dG'].append(mm_dG)
            results['mmgbsa_err'].append(mm_dG - exp_dG)

    df_dir = os.path.join('analysis', setname)
    if not os.path.isdir(df_dir): os.makedirs(df_dir)

    df = pd.DataFrame(data=results, dtype='object')
    df.to_csv(f'{df_dir}/results.csv', index=False)

    df_exp = df.sort_values('exp_dG', inplace=False)
    df_exp.to_csv(f'{df_dir}/results_sorted_exp.csv', index=False)

    df_dock = df.sort_values('dock_dG', inplace=False)
    df_dock.to_csv(f'{df_dir}/results_sorted_dock.csv', index=False)

    df_apnet = df.sort_values('apnet_dG', inplace=False)
    df_apnet.to_csv(f'{df_dir}/results_sorted_apnet.csv', index=False)

    if do_mmgbsa:
        df_mmgbsa = df.sort_values('mmgbsa_dG', inplace=False)
        df_mmgbsa.to_csv(f'{df_dir}/results_sorted_mmgbsa.csv', index=False)

    summary = open(f'{df_dir}/summary.txt', 'w')
    summary.write(f'\tSUMMARY OF RESULTS FOR SYSTEM {setname.upper()}\n\n')

    dock_err = results['dock_err']
    write_summary('PRIME', dock_err, results['dock_dG'], results['exp_dG'], summary)

    apnet_err = results['apnet_err']
    write_summary('AP-NET-DG', apnet_err, results['apnet_dG'], results['exp_dG'], summary)

    if do_mmgbsa:
        mmgbsa_err = results['mmgbsa_err']
        write_summary('MMGBSA', mmgbsa_err, results['mmgbsa_dG'], results['exp_dG'], summary)