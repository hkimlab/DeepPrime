import os, sys, time, warnings
import multiprocessing as mp
import numpy as np

warnings.filterwarnings('ignore')

from src.biofeat import *
from src.dspcas9 import calculate_DeepSpCas9_score
from src.dprime import calculate_deepprime_score

np.set_printoptions(threshold=sys.maxsize)

time_start = time.time()

## system config ##
sBASE_DIR   = os.getcwd()
nCore_max   = os.cpu_count()
CORE_CNT    = 4

def main_run(analysistag):
    work_dir    = '%s/data/%s' % (sBASE_DIR, analysistag)
    seqs        = [readline.strip('\n').split(',') for readline in open('%s/seqs.txt' % work_dir, 'r')]
    options     = [readline.strip('\n').split(',') for readline in open('%s/options.txt' % work_dir, 'r')][0]

    pe_type     = options[0]
    pbs_min     = int(options[1])
    pbs_max     = int(options[2])
    rtt_max     = int(options[3])

    ## Parameters
    dict_params = {'nAltIndex': 60,                  # 60nts --- Alt --- 60nts *0-based
                   'bTest': 0,                       # 0 or 1 / True or False
                   'PBS_range': [pbs_min, pbs_max],  # Range limit = 1 - 17
                   'RTT_max':   rtt_max,               # Range limit = 40
                   'PE_system': pe_type              # PE2 / NRCH_PE2 / PE2max
                   }

    print('\n\nStart: %s - %s\n\n' % (analysistag, dict_params['PE_system']))

    ## Load input file
    format_input(work_dir, seqs, options)

    df_input   = pd.read_csv('%s/input.csv' % work_dir)
    n_input    = len(df_input)
    bins       = [[int(n_input * (i + 0) / CORE_CNT), int(n_input * (i + 1) / CORE_CNT)] for i in range(CORE_CNT)]

    parameters = []
    for start, end in bins:
        df_subsets = df_input[start:end]
        parameters.append([analysistag, work_dir, df_subsets, dict_params])
        print('Input subset range:', start, end)
    # loop END: nStart, nEnd


    ## Multiprocessing
    p = mp.Pool(CORE_CNT)
    p.map_async(mp_processor, parameters).get()
    p.close()
    p.join()

    print('\n--- All multiprocessing finished ---\n')
    print('%.3f sec\n\n' % (time.time() - time_start))
# def END: main_run

def format_input(work_dir, seqs, options, target='on'):

    editinfo    = options[4]

    list_output = ['pegRNA%s,%s,%s,%s' % (i + 1, seqs[i][0], seqs[i][1], editinfo) for i in range(len(seqs))]

    out         = '%s/input%s.csv' % (work_dir, target if target == 'off' else '')
    outfile     = open(out, 'w')
    outfile.write('#id,%sseq,EDseq,EditInfo\n' % target)
    for out in list_output:
        outfile.write('%s\n' % out)
    outfile.close()
# def END: format_input


def off_run (sANALYSISTAG):

    work_dir      = '%s/data/%s' % (sBASE_DIR, sANALYSISTAG)
    options       = [readline.strip('\n').split(',') for readline in open('%s/options.txt' % work_dir, 'r')][0]
    list_pegRNAs  = [readline.strip('\n').split(',') for readline in open('%s/pegRNA.txt' % work_dir, 'r')][0]
    list_offseqs  = [readline.strip('\n').split(',') for readline in open('%s/offseq.txt' % work_dir, 'r')][0]

    ## Step 1: make input ##
    make_input      (work_dir, list_pegRNAs, list_offseqs)

    ## Step 2: offtarget test ##
    off_target_test (work_dir)
#def END: off_run


def make_input (work_dir, list_pegRNAs, list_offseqs):
        dict_target, h = get_csv_data (work_dir, list_pegRNAs)
        out            = '%s/input_off.csv' % work_dir
        outfile        = open(out, 'w')
        outfile.write('ID,WT74_ref,%s\n' % ','.join(h))

        for pegID in dict_target:
            for i, offseq in enumerate(list_offseqs):
                outfile.write('%s-off%s,%s,%s\n' % (pegID, (i+1) , offseq, ','.join(dict_target[pegID])))
        outfile.close()
#def END: make_input


def off_target_test (work_dir):

    df = pd.read_csv('%s/input_off.csv' % work_dir)

    df['Offtarget_Score'] = calculate_deepprime_score(sBASE_DIR, df, 'PE2_offtarget')
    sorted_df = df.sort_values(by=['Offtarget_Score'], ascending=False)
    print(sorted_df)
    sorted_df.to_csv('%s/offtarget_output.csv' % (work_dir), index=False)
#def END: off_target_test


def get_csv_data (work_dir, list_pegRNAs):

    dict_output = {}
    output_dir  = '%s/output/results' % work_dir

    for pegID in list_pegRNAs:

        if pegID not in dict_output:
            dict_output[pegID] = ''

        pegRNA, ID = pegID.split('.')
        in_csv     = '%s/%s.csv' % (output_dir, pegRNA)

        dict_info  = {readline.strip('\n').split(',')[0]:readline.strip('\n').split(',')[1:] for readline in open(in_csv)}
        header     = dict_info['ID']
        dict_output[pegID] = dict_info[pegID]
    #loop END: pegID

    return dict_output, header
#def END: get_csv_data



def mp_processor(parameters):
    analysistag, work_dir, df_subsets, dict_params = parameters

    out_dir  = '%s/output' % work_dir
    os.makedirs(out_dir, exist_ok=True)

    for idx in df_subsets.index:
        df_temp          = df_subsets.loc[idx]
        id               = df_temp[0]
        wtseq            = df_temp[1].upper()
        editedseq        = df_temp[2].upper()
        edit             = df_temp[3]
        deep_prime(analysistag, id, wtseq, editedseq, edit, dict_params)

    # loop END: idx

    print('Done')
# def END: mp_processor


def deep_prime(analysistag, id, wtseq, editedseq, edit, dict_params):

    nAltIndex   = dict_params['nAltIndex']
    bTest       = dict_params['bTest']
    pbs_range   = dict_params['PBS_range']
    rtt_max     = dict_params['RTT_max']
    pe_system   = dict_params['PE_system']
    edit_type   = edit[:-1]
    edit_len    = int(edit[-1])

    output      = '%s/data/%s/output' % (sBASE_DIR, analysistag)
    os.makedirs(output, exist_ok=True)

    ## FeatureExtraction Class
    cFeat = FeatureExtraction()

    cFeat.input_id = id
    cFeat.get_input(wtseq, editedseq, edit_type, edit_len)

    cFeat.get_sAltNotation(nAltIndex)
    cFeat.get_all_RT_PBS(nAltIndex, nMinPBS=pbs_range[0] - 1, nMaxPBS=pbs_range[1], nMaxRT=rtt_max, pe_system=pe_system)
    cFeat.make_rt_pbs_combinations()
    cFeat.determine_seqs()
    cFeat.determine_secondary_structure()

    df = cFeat.make_output_df(bTest)

    if len(df) == 0:  # Empty DataFrame = No PAM found
        return [cFeat.input_id, 0, 0, 0, 0, 0, 0, 0], 'NULL'

    else:
        list_Guide30 = [WT74[:30] for WT74 in df['WT74_On']]

        df['DeepSpCas9_score'] = calculate_DeepSpCas9_score(sBASE_DIR, list_Guide30)
        df['DeepPrime_score']  = calculate_deepprime_score(sBASE_DIR, df, pe_system)
        sorted_df = df.sort_values(by=['DeepPrime_score'], ascending=False)

        df_top_pegs = sorted_df[:10]

        ## Save result file
        # df.to_feather('%s/%s.feather' % (outdir, sID))
        sorted_df.to_csv('%s/%s.csv' % (output, id), index=False)

        dp_score        = df.DeepPrime_score
        tot_pegRNAs     = len(dp_score)
        ave_DPScore     = np.mean(dp_score)
        ave_DPStop5     = np.mean(dp_score.sort_values(ascending=False).head(5))
        list_stat       = [cFeat.input_id, tot_pegRNAs, ave_DPScore, ave_DPStop5]

        dict_top        = df_top_pegs.to_dict('index')
        dict_topdesign  = get_oligos(dict_top)
        outf            = open('%s/%s.top_designs.csv' % (output, id), 'w')
        header          = '#id,wtseq, edseq, spacer_top, spacer_bot, ext_top, ext_bot, dpcas, dpprime\n'
        outf.write(header)
        for id in dict_topdesign:
            outf.write('%s,%s\n' % (id, ','.join([str(data) for data in dict_topdesign[id]])))
        #loop END: id
        outf.close()
        return list_stat, df_top_pegs
    #if END:
# def END: deep_prime

def get_oligos (dict_top):

    dict_output = {}

    for index in dict_top:
        id          = dict_top[index]['ID']
        wtseq       = dict_top[index]['WT74_On']
        edseq       = dict_top[index]['Edited74_On']
        gN19        = dict_top[index]['gN19']
        RTPBS       = dict_top[index]['Edited74_On'].strip('x')
        dpcas       = dict_top[index]['DeepSpCas9_score']
        dpprime     = dict_top[index]['DeepPrime_score']

        spacer_top  = 'cacc%sgtttt' % gN19
        spacer_bot  = 'ctctaaaac%s' % reverse_complement(gN19)
        ext_top     = 'gtgc%s' % RTPBS
        ext_bot     = 'aaaa%s' % reverse_complement(RTPBS)

        if id not in dict_output:
            dict_output[id] = ''
        dict_output[id] = [wtseq, edseq, spacer_top, spacer_bot, ext_top, ext_bot, dpcas, dpprime]
    #loop END: index

    return dict_output
# def END: get_oligos


def main():
    print('Usage')
    print('main.py main_run <AnalysisTag>')


# def END: main


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        function_name = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys():
            locals()[function_name](*function_parameters)
        else:
            sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    # if END: len(sys.argv)
# if END: __name__
