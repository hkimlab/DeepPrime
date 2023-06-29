from src.biofeat import *
from src.utils import print_logo
print_logo()

import os, sys, time, argparse, glob
import numpy as np
import multiprocessing as mp
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

from src.dspcas9 import calculate_DeepSpCas9_score
from src.dprime import calculate_deepprime_score


np.set_printoptions(threshold=sys.maxsize)

def deepprime_pipeline(input_file, sample_name, pe_type, cell_type, pbs_min, pbs_max, jobs, progress):

    time_start = time.time()

    ## system config
    nCore_max = os.cpu_count()
    nCores    = jobs

    ## Parameters
    dict_params = { 'nAltIndex' : 60,                  # 60nts --- Alt --- 60nts *0-based
                    'bTest'     : 0,                   # 0 or 1 / True or False
                    'PBS_range' : [pbs_min, pbs_max],  # Range limit = 1 - 17
                    'RTT_max'   : 40,                  # Range limit = 40
                    'PE_system' : pe_type,             # PE2 / NRCH_PE2 / PE2max
                    'Cell_type' : cell_type,
                    }
    
    print('[Info] Start DeepPrime pipeline\nID: %s\nPE system: %s\n' % (sample_name, dict_params['PE_system']))

    ## Load input file
    df_input = pd.read_csv(input_file)
    n_input  = len(df_input)
    list_nBins = [[int(n_input * (i + 0) / nCores), int(n_input * (i + 1) / nCores)] for i in range(nCores)]
    list_sParameters = []

    # Make output directory
    sBase_DIR=os.getcwd()
    out_dir = '%s/%s_output/%s' % (sBase_DIR, sample_name, dict_params['PE_system'])
    os.makedirs(out_dir, exist_ok=True)

    for nStart, nEnd in list_nBins:
        df_sSubSplits = df_input[nStart:nEnd]
        list_sParameters.append([sample_name, df_sSubSplits, dict_params, nStart, nEnd, out_dir, progress])
        print('Input subset range:', nStart, nEnd)
    
    print('\n[Info] Total number of processing cores: %s' % jobs)

    ## Multiprocessing
    p = mp.Pool(nCores)
    p.map_async(worker, list_sParameters).get()
    p.close()
    p.join()

    # Make summary files
    df_stat = make_merged_file(out_dir, '*stat.csv', remove_temp=True)
    df_stat.to_csv('%s/Statistics.csv' % out_dir, index=False)

    df_top4 = make_merged_file(out_dir, 'Top4_pegRNAs_*.csv', remove_temp=True)
    df_top4.to_csv('%s/Top4_pegRNAs.csv' % out_dir, index=False)

    print('\n--- All multiprocessing finished ---\n')
    print('%.3f sec\n\n' % (time.time() - time_start))

# def END: main


def worker(list_sParameters: list):
    sample_name, df_sSubSplits, dict_params, nStart, nEnd, out_dir, progress = list_sParameters

    list_output = []
    df_top4_lib = pd.DataFrame()

    for idx in df_sSubSplits.index:
        df_temp = df_sSubSplits.loc[idx]
        
        sID      = df_temp[0]
        Ref_seq  = df_temp[1].upper()
        ED_seq   = df_temp[2].upper()
        sAlt     = df_temp[3]
        
        # DeepPrime-Off model
        if dict_params['PE_system'] == 'PE-off':
            print('\nDeepPrime-Off model is not available due to pipeline renewal')
            print('Will be available soon. Please contact gsyu93@gmail.com if you need help.')
            sys.exit()
        
        else: result, DPStop4 = deep_prime(sample_name, sID, Ref_seq, ED_seq, sAlt, out_dir, dict_params)
        
        list_output.append(result)
        
        if progress: print('Processing: %s' % sID)

        # progress = round(100 * (idx / len(df_sSubSplits)), 1)
        # print('Processing: %d%% - %s' % (progress, sID))

        df_top4_lib = pd.concat([df_top4_lib, DPStop4], ignore_index=True)

    columnes = ['ID', 'Total_pegRNAs', 'Average_DP_score', 'Best_pegRNA_score', 'Average_Top4_DP_score',
                'Over30_pegRNAs', 'Over20_pegRNAs', 'Over10_pegRNAs', 'Over5_pegRNAs']

    df_stat = pd.DataFrame(list_output, columns=columnes)
    df_stat.to_csv('%s/%s_%s_%d_%d-%dstat.csv' % (out_dir, sample_name, dict_params['PE_system'], idx, nStart, nEnd), index=False)
    df_top4_lib.to_csv('%s/Top4_pegRNAs_%s_%s_%d-%d.csv' % (out_dir, sample_name, dict_params['PE_system'], nStart, nEnd), index=False)

# def END: mp_processor

def deep_prime(sample_name: str,
                sID: str, 
                Ref_seq: str, 
                ED_seq: str, 
                sAlt: str, 
                output_dir: str,
                dict_params=None, 
                ):

    ## Default parameters
    default_params = { 'nAltIndex'   : 60,  # 60nts --- Alt --- 60nts *0-based
                    'bTest'       : 0,
                    'PBS_range'   : [1, 17],
                    'RTT_max'     : 40,
                    'PE_system'   : 'PE2'
                    }

    if dict_params: parameters = dict_params
    else:           parameters = default_params

    nAltIndex   = parameters['nAltIndex']
    bTest       = parameters['bTest']
    pbs_range   = parameters['PBS_range']
    rtt_max     = parameters['RTT_max']
    pe_system   = parameters['PE_system']
    cell_type   = parameters['Cell_type']

    edit_type   = sAlt[:-1]
    edit_len    = int(sAlt[-1])

    results_dir = '%s/results' % output_dir
    os.makedirs(results_dir, exist_ok=True)

    ## FeatureExtraction Class
    cFeat = FeatureExtraction()

    cFeat.input_id = sID
    cFeat.get_input(Ref_seq, ED_seq, edit_type, edit_len)

    cFeat.get_sAltNotation(nAltIndex)
    cFeat.get_all_RT_PBS(nAltIndex, nMinPBS=pbs_range[0]-1, nMaxPBS=pbs_range[1], nMaxRT=rtt_max, pe_system=pe_system)
    cFeat.make_rt_pbs_combinations()
    cFeat.determine_seqs()
    cFeat.determine_secondary_structure()

    df = cFeat.make_output_df(bTest)

    if len(df) == 0: # Empty DataFrame = No PAM found
        return [cFeat.input_id, 0, 0, 0, 0, 0, 0, 0, 0], pd.DataFrame()

    else:
        list_Guide30 = [WT74[:30] for WT74 in df['WT74_On']]
        df['DeepSpCas9_score'] = calculate_DeepSpCas9_score(os.getcwd(), list_Guide30)
        df['DeepPrime_score']  = calculate_deepprime_score(df, pe_system, cell_type)

        ## Save result file
        df.to_csv('%s/%s.csv' % (results_dir, sID), index=False)

        dp_score = df.DeepPrime_score

        tot_pegRNAs = len(dp_score)
        ave_DPScore = np.mean(dp_score)
        ave_DPStop1 = np.mean(dp_score.sort_values(ascending=False).head(1))
        ave_DPStop4 = np.mean(dp_score.sort_values(ascending=False).head(4))

        DPStop4 = df.sort_values(by='DeepPrime_score', ascending=False).head(4)

        over_5  = dp_score[dp_score >= 5]
        over_10 = over_5[over_5 >= 10]
        over_20 = over_10[over_10 >= 20]
        over_30 = over_20[over_20 >= 30]

        list_stat = [cFeat.input_id, tot_pegRNAs, ave_DPScore, ave_DPStop1, ave_DPStop4, len(over_30), len(over_20), len(over_10), len(over_5)]

        return list_stat, DPStop4

# def END: deep_prime

def make_merged_file(file_path:str, name_pattern:str, remove_temp=True):
    files = glob.glob(os.path.join(file_path + '/', name_pattern))

    list_output_df = []

    for file in files:
        mer_out = pd.read_csv(file)
        list_output_df.append(mer_out)

        if remove_temp: os.remove(file)

    df_merged = pd.concat(list_output_df, ignore_index=True)

    return df_merged


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input_file", type=str, 
                        help="Input file containing target sequence and edit type")
    
    parser.add_argument("-n", "--name", type=str, 
                        help="Sample name for your input (default='Sample')", default='Sample')
    
    parser.add_argument("-p", "--pe_type", type=str, 
                        choices=['PE2', 'PE2max', 'PE2max-e', 'PE4max', 'PE4max-e', 'NRCH_PE2', 'NRCH_PE2max', 'NRCH_PE4max', 'PE-off'], 
                        help="PE type parameter (default=PE2max)", default='PE2max')
    
    parser.add_argument("--cell_type", 
                        type=str, 
                        choices=['HEK293T', 'A549', 'DLD1', 'HCT116', 'HeLa', 'MDA-MB-231', 'NIH3T3'],
                        help="Cell type parameter. (default=HEK293T)", default='HEK293T')
    
    parser.add_argument("--pbs_min", type=int, 
                        help="PBS minimum length parameter (default=1)", default=1)
    
    parser.add_argument("--pbs_max", type=int, 
                        help="PBS maximun length parameter (default=17)", default=17)
    
    parser.add_argument("--jobs", type=int, 
                        help="Number of cores for computing (default=1)", default=1)
    
    parser.add_argument("--progress", action='store_true',
                        help="Show processing message")

    # 입력된 파라미터를 받아와 사용합니다.
    args = parser.parse_args()

    input_file  = args.input_file
    sample_name = args.name
    pe_type     = args.pe_type
    cell_type   = args.cell_type
    pbs_min     = args.pbs_min
    pbs_max     = args.pbs_max
    jobs        = args.jobs
    progress    = args.progress

    if pe_type == 'PE-off':
        print(
'''
#####################################################################
[Notice] DeepPrime-Off is not available at this version of pipeline.

DeepPrime-Off model is not available due to pipeline renewal.
You can use model check point file (models/DeepPrime/DeepPrime_off).

The pipeline version Will be available soon.
Please contact gsyu93@gmail.com, if you need help.

Thank you.
#####################################################################
'''
              )
        
        sys.exit()

    deepprime_pipeline(input_file, sample_name, pe_type, cell_type, pbs_min, pbs_max, jobs, progress)

# if END: __name__
