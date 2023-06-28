import os, sys, time, warnings
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator
from scipy import stats
from sklearn.preprocessing import minmax_scale


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
        deep_prime(analysistag, id, wtseq, editedseq, edit, dict_params, out_dir)

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
    top         = 5

    results    = '%s/results' % outdir
    plots      = '%s/matplot' % outdir
    os.makedirs(results, exist_ok=True)
    os.makedirs(plots, exist_ok=True)

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
        df['DeepPrime_score'], df['Zscores']  = calculate_deepprime_score(sBASE_DIR, df, pe_system)

        sorted_df = df.sort_values(by=['DeepPrime_score'], ascending=False)

        df_top_pegs = sorted_df[:top]
        dict_top    = df_top_pegs.to_dict('index')


        list_scores     = sorted(df['Zscores'])
        list_minmax     = minmax_scale(list_scores)

        list_topdesign  = get_oligos(plots, list_scores, dict_top)

        rank_plots (plots, list_scores, list_minmax, list_topdesign)

        header          = 'id,guide,pbsrtt,wtseq,edseq,spacer_top,spacer_bot,ext_top,ext_bot,dpcas,dpprime,pbslen,rttlen,edittype,zscore,percentile,plotpath'
        outf = open('%s/%s.top_designs.csv' % (results, id), 'w')
        outf.write('%s\n' % header)
        for data in list_topdesign:
            out = ','.join([str(col) for col in data])
            outf.write('%s\n' % out)
        #loop END: id
        outf.close()

        ## Save result file
        # df.to_feather('%s/%s.feather' % (outdir, sID))
        sorted_df.to_csv('%s/%s.csv' % (results, id), index=False)
        return list_topdesign
    #if END:
# def END: deep_prime

def get_oligos (plotphp, list_scores, dict_top):

    list_output = []
    for index in dict_top:
        id          = dict_top[index]['ID']
        wtseq       = dict_top[index]['WT74_On']
        edseq       = dict_top[index]['EditedwNote']
        gN19        = dict_top[index]['gN19']
        RTPBS       = dict_top[index]['Edited74_On'].strip('x')
        dpcas       = dict_top[index]['DeepSpCas9_score']
        dpprime     = dict_top[index]['DeepPrime_score']
        pbslen      = dict_top[index]['PBSlen']
        rttlen      = dict_top[index]['RTlen']
        edit        = dict_top[index]['AltKey']
        zscore      = dict_top[index]['Zscores']
        percentile  = '%sth' % (int(stats.percentileofscore(list_scores, zscore)))

        plotpath    = '%s/%s.png' % (plotphp, id)

        spacer_top  = 'cacc%sgtttt' % gN19
        spacer_bot  = 'ctctaaaac%s' % reverse_complement(gN19)
        ext_top     = 'gtgc%s' % reverse_complement(RTPBS)
        ext_bot     = 'aaaa%s' % RTPBS

        list_output.append([id, gN19, RTPBS, wtseq, edseq, spacer_top, spacer_bot, ext_top, ext_bot, dpcas, dpprime, pbslen, rttlen, edit, zscore, percentile, plotpath])
    #loop END: index

    return list_output
# def END: get_oligos


def rank_plots (plots, list_scores, list_percents, list_top):
    for data in list_top:
        pegID      = data[0]
        zscore     = data[14]
        rank       = sorted(list_scores + [zscore]).index(zscore) # used as index for line
        outplot    = '%s/%s.png' % (plots, pegID)

        print('MATPLOTLIB - Rank Plot - %s' % pegID)
        matplot_scatter (list_scores, list_percents, outplot, zscore)
    #loop END: data
# def END: rank_plots


def matplot_scatter(list_scores, list_percents, outf, rank):
    list_X, list_Y  = [], []
    for x, y in zip(list_scores, list_percents):
        list_X.append(x)
        list_Y.append(y*100)

    assert len(list_X) == len(list_Y)
    ### Figure Size ###
    FigWidth         = 12
    FigHeight        = 5

    OutFig  = plt.figure(figsize=(FigWidth, FigHeight))
    SubPlot = OutFig.add_subplot(111)

    ### Marker ###########
    Red             = 0
    Green           = 0
    Blue            = 0
    MarkerSize      = 25
    Circle          = 'o'
    DownTriangle    = 'v'
    #######################

    ### Log Start Point ###
    LimitThresh     = 10
    #######################

    ### Axis Range ########
    Xmin = min(list_scores)
    Xmax = max(list_scores)
    Ymin = 0
    Ymax = 100
    ########################

    ### Tick Marks #########
    TickScale     = 1
    MajorTickSize = 10
    MinorTickSize = 5

    plt.xlim(xmin=Xmin, xmax=Xmax)
    plt.ylim(ymin=Ymin, ymax=Ymax)
    #plt.xscale('symlog', linthreshx=LimitThresh)
    #plt.yscale('symlog', linthreshy=LimitThresh)

    #plt.axes().xaxis.set_minor_locator(MinorSymLogLocator(TickScale))
    #plt.axes().yaxis.set_minor_locator(MinorSymLogLocator(TickScale))

    plt.tick_params(which='major', length=MajorTickSize)
    plt.tick_params(which='minor', length=MinorTickSize)

    #Set xticks as percentiles
    #p = np.array([0.0, 25.0, 50.0, 75.0, 100.0])
    #plt.yticks((len(list_Y) - 1) * p / 100., map(str, p))

    plt.vlines(rank, ymin=Ymin, ymax=Ymax, colors='red')
    plt.xlabel('DeepPrime Score')
    plt.ylabel('Score Rank')
    SubPlot.scatter(list_X, list_Y, marker=Circle, color='black', s=MarkerSize)
    OutFig.savefig(outf)

#def END: matplot



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
