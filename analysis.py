import os, sys, argparse
import pandas as pd
import numpy as np
import schrodinger.structure as structure

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

def write_summary(name, errs, file):
    me = np.mean(errs)
    rmse = np.sqrt(np.mean(np.square(errs)))
    mae = np.mean(np.abs(errs))

    file.write(f'\t{name.upper()} ERRORS => ME: {me}, RMSE: {rmse}, MAE: {mae}\n\n')

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
        'exp_dG' : [],
        'prime_dG' : [],
        'prime_err' : [],
        'apnet_dG' : [],
        'apnet_err' : [],
    }

    if do_mmgbsa:
        results['mmgbsa_dG'] = []
        results['mmgbsa_err'] = []

    dg_data = pd.read_pickle(f'{dg_path}/datasets/{setname}/dimers.pkl')
    for n in range(len(dg_data)):
        results['system'].append(dg_data['system'][n])
        results['exp_dG'].append(to_dG(dg_data['label'][n]))

    prime_datapath = f'docking/{setname}/dockjob_pv.maegz'
    for n, st in enumerate(structure.StructureReader(prime_datapath)):
        if n == 0: continue
        prime_dG = st.property['r_i_glide_gscore']
        ref_dG = results['exp_dG'][n-1]
        results['prime_dG'].append(prime_dG)
        results['prime_err'].append(prime_dG - ref_dG)

    dg_pred = pd.read_csv(f'apnetdg/{setname}/preds.csv')
    for n in range(len(dg_pred)):
        pred_dG = to_dG(dg_pred['AP_net_dG_pKi'][n])
        ref_dG = results['exp_dG'][n]
        results['apnet_dG'].append(pred_dG)
        results['apnet_err'].append(pred_dG - ref_dG)

    if do_mmgbsa:
        mm_res = pd.read_csv(f'mmgbsa/{setname}/output.csv')
        for n in range(len(mm_res)):
            mm_dG = mm_res['r_psp_MMGBSA_dG_Bind'][n]
            ref_dG = results['exp_dG'][n]
            results['mmgbsa_dG'].append(mm_dG)
            results['mmgbsa_err'].append(mm_dG - ref_dG)

    df_dir = os.path.join('analysis', setname)
    if not os.path.isdir(df_dir): os.makedirs(df_dir)

    df = pd.DataFrame(data=results, dtype='object')
    df.to_csv(f'{df_dir}/results.csv', index=False)

    df_exp = df.sort_values('exp_dG', inplace=False)
    df_exp.to_csv(f'{df_dir}/results_sorted_exp.csv', index=False)

    df_prime = df.sort_values('prime_dG', inplace=False)
    df_prime.to_csv(f'{df_dir}/results_sorted_prime.csv', index=False)

    df_apnet = df.sort_values('apnet_dG', inplace=False)
    df_apnet.to_csv(f'{df_dir}/results_sorted_apnet.csv', index=False)

    if do_mmgbsa:
        df_mmgbsa = df.sort_values('mmgbsa_dG', inplace=False)
        df_mmgbsa.to_csv(f'{df_dir}/results_sorted_mmgbsa.csv', index=False)

    summary = open(f'{df_dir}/summary.txt', 'w')
    summary.write(f'\tSUMMARY OF RESULTS FOR SYSTEM {setname.upper()}\n\n')

    prime_err = results['prime_err']
    write_summary('PRIME', prime_err, summary)

    apnet_err = results['apnet_err']
    write_summary('AP-NET-DG', apnet_err, summary)

    if do_mmgbsa:
        mmgbsa_err = results['mmgbsa_err']
        write_summary('MMGBSA', mmgbsa_err, summary)