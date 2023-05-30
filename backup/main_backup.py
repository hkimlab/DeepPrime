import os, sys, time, warnings
from scipy import stats
from sklearn.preprocessing import minmax_scale

import numpy as np
import multiprocessing as mp
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator

from Bio import Entrez, GenBank, SeqIO
Entrez.email = "gsyu93@gmail.com"
import xml.etree.ElementTree as ET

warnings.filterwarnings('ignore')

from src.biofeat import *
from src.dspcas9 import calculate_DeepSpCas9_score
from src.dprime import calculate_deepprime_score

np.set_printoptions(threshold=sys.maxsize)

time_start = time.time()

## system config ##
sBASE_DIR = os.getcwd()
nCore_max = os.cpu_count()
CORES     = 1
ALT_INDEX = 60  # 60nts --- Alt --- 60nts *0-based
bTEST     = 0   # 0 or 1

def main_run(analysistag):

    ## Set Parameters ##
    work_dir    = '%s/data/%s' % (sBASE_DIR, analysistag)
    wtseqs      = [readline.strip('\n').replace(' ', '').split(',') for readline in open('%s/wtseq.txt' % work_dir, 'r')][0]
    editedseqs  = [readline.strip('\n').replace(' ', '').split(',') for readline in open('%s/editedseq.txt' % work_dir, 'r')][0]
    options     = [readline.strip('\n').split(',') for readline in open('%s/options.txt' % work_dir, 'r')][0]
    pe_type     = options[0]
    pbs_min     = int(options[1])
    pbs_max     = int(options[2])
    rtt_max     = int(options[3])
    editinfo    = options[4]

    dict_params = {
                   'PBS_range': [pbs_min, pbs_max],  # Range limit = 1 - 17
                   'RTT_max':   rtt_max,             # Range limit = 40
                   'PE_system': pe_type              # PE2 / NRCH_PE2 / PE2max
                   }
    print('\n\nStart: %s - %s\n\n' % (analysistag, dict_params['PE_system']))

    ## Load input file
    format_input(work_dir, wtseqs, editedseqs, editinfo)
    df_input   = pd.read_csv('%s/input.csv' % work_dir)
    n_input    = len(df_input)
    list_nBins = [[int(n_input * (i + 0) / CORES), int(n_input * (i + 1) / CORES)] for i in range(CORES)]

    parameters = []
    for nStart, nEnd in list_nBins:
        df_sSubSplits = df_input[nStart:nEnd]
        parameters.append([analysistag, work_dir, df_sSubSplits, dict_params])
        print('Input subset range:', nStart, nEnd)
    # loop END: nStart, nEnd

    ## Multiprocessing
    p = mp.Pool(CORES)
    p.map_async(mp_processor, parameters).get()
    p.close()
    p.join()

    print('\n--- All multiprocessing finished ---\n')
    print('%.3f sec\n\n' % (time.time() - time_start))
# def END: main_run


def check_clinvar (analysistag):
    dict_sCLNVCKey = {'Insertion': 'insertion',
                      'Duplication': 'insertion',

                      'Deletion': 'deletion',

                      'Inversion': 'substitution',
                      'single_nucleotide_variant': 'substitution',
                      'single nucleotide variant': 'substitution'}

    genome          = 'GRCh38'
    work_dir        = '%s/data/%s' % (sBASE_DIR, analysistag)
    clivvarfile     = '%s/clinvar_run.txt' % work_dir

    seqinfo         = get_seqs(genome, clivvarfile, ALT_INDEX, dict_sCLNVCKey)

    try:
        if seqinfo.startswith('Invalid'):
            print(seqinfo)
            return

    except AttributeError:
        wtseq    = seqinfo[0]
        edseq    = seqinfo[1]
        editinfo = seqinfo[2]

        outf = open('%s/wtseq.txt' % work_dir, 'w')
        # outf.write('%s\n' % ','.join([seq for seq in wtseqs]))
        outf.write('%s\n' % wtseq)
        outf.close()

        outf = open('%s/editedseq.txt' % work_dir, 'w')
        # outf.write('%s\n' % ','.join([seq for seq in edseqs]))
        outf.write('%s\n' % edseq)
        outf.close()

        options = [readline.strip('\n').split(',') for readline in open('%s/options.txt' % work_dir, 'r')][0]
        updated = options[:-1] + [editinfo]
        outf    = open('%s/options.txt' % work_dir, 'w')
        outf.write('%s\n' % ','.join([para for para in updated]))
        outf.close()

        outf = open('%s/check.txt' % work_dir, 'w')
        outf.write('%s\n' % (','.join(seqinfo)))
        outf.close()

        temp_dir    = '%s/temp' % sBASE_DIR
        os.makedirs(temp_dir, exist_ok=True)

        svgfile     = '%s/clinvarcheck.svg' % work_dir
        pngfile     = '%s/clinvarcheck.png' % work_dir
        #os.system('cp %s %s' % (temppng, clinvarplot))
        #svg_clinvar (svgfile, pngfile, wtseq, edseq, editinfo)

        print('PASS')
    #try END:
#def END: check_runtype


def get_seqs(genome, clivvarfile, altindex, dict_sCLNVCKey):

    inputID, run = [readline.strip('\n').replace(' ', '').split(',') for readline in open(clivvarfile)][0]

    if inputID.startswith('VCV'):
        query  = inputID.split('.')[0]
        handle = Entrez.efetch(db='clinvar', id=query, rettype='vcv')
    else:
        handle = Entrez.efetch(db='clinvar', id=inputID, rettype='vcv', is_varationid='true', from_esearch="true")
    #if END:

    result    = ET.parse(handle)
    handle.close()
    root      = result.getroot()

    var       = root.findall('./VariationArchive')[0]
    vartype   = var.attrib['VariationType']
    var_loc   = root.findall('./VariationArchive/InterpretedRecord/SimpleAllele/Location/SequenceLocation')
    list_info = []

    for info in var_loc:
        if info.attrib['Assembly'] != genome: continue
        list_info.append([info.attrib['Accession'], int(info.attrib['start']), info.attrib['referenceAlleleVCF'], info.attrib['alternateAlleleVCF'], info.attrib['variantLength']])
    #loop END: info
    assert len(list_info) == 1

    accID, start, ref_nt, alt_nt, altlen = list_info[0]
    chr_seq_fetch = Entrez.efetch(db="nucleotide",
                                  id=accID,
                                  rettype="fasta",
                                  strand=1,
                                  seq_start=start - altindex,
                                  seq_stop=start + altindex)

    record           = SeqIO.read(chr_seq_fetch, "fasta")
    chr_seq_fetch.close()

    editinfo         = determine_alttype_altlen(vartype, ref_nt, alt_nt, dict_sCLNVCKey)
    if editinfo.startswith('Invalid'): return editinfo

    seq              = str(record.seq)

    if run == 'model':
        wtseq  = seq
        edseq  = seq[:nAltIndex] + alt_nt + seq[nAltIndex + len(ref_nt):]

    else:
        wtseq  = seq[:nAltIndex] + alt_nt + seq[nAltIndex + len(ref_nt):]
        edseq  = seq
    #if END:

    return wtseq, edseq, editinfo
#def END: get_seqs


def determine_alttype_altlen (vartype, ref_nt, alt_nt, dict_sCLNVCKey):

    if vartype in ['Microsatellite', 'Indel']:
        if len(ref_nt) > len(alt_nt):
            alttype = 'deletion'
            altlen  = len(ref_nt) - 1

            if ref_nt[0] != alt_nt[0]:
                return 'Invalid Ref-Alt nucleotides:\\nFound: %s -> %s' % (ref_nt, alt_nt) # ex) GAGA -> TCT

        elif len(ref_nt) < len(alt_nt):
            alttype = 'insertion'
            altlen  = len(alt_nt) - 1

            if ref_nt[0] != alt_nt[0]:
                return 'Invalid Ref-Alt nucleotides:\\nFound: %s -> %s' % (ref_nt, alt_nt) # ex) GAGA -> TCT
        else:
            alttype = 'substitution'
            altlen  = len(alt_nt)
    else:
        try: alttype = dict_sCLNVCKey[vartype]
        except KeyError: return 'Invalid Variant Type:\\nFound: %s' % (vartype)

        if alttype == 'deletion':
            altlen = len(ref_nt) - 1

        elif alttype == 'insertion':
            altlen = len(alt_nt) - 1
        else:
            altlen = len(alt_nt)

    if altlen > 3: return 'Invalid Variant Size:\\nFound: %snt' % (altlen)

    editinfo = '%s%s' % (alttype[:3], altlen)

    return editinfo
#def END: determine_alttype_altlen


def format_input(work_dir, wtseqs, editedseqs, editinfo, target='on'):

    list_output = ['pegRNA%s,%s,%s,%s' % (i + 1, wtseqs[i], editedseqs[i], editinfo) for i in range(len(wtseqs))]
    out         = '%s/input%s.csv' % (work_dir, target if target == 'off' else '')
    outfile     = open(out, 'w')
    outfile.write('#id,%sseq,EDseq,EditInfo\n' % target)
    for out in list_output:
        outfile.write('%s\n' % out)
    outfile.close()
# def END: format_input


def off_run (analysistag):

    work_dir      = '%s/data/%s' % (sBASE_DIR, analysistag)
    list_pegRNAs  = [readline.strip('\n').split(',') for readline in open('%s/pegRNA.txt' % work_dir, 'r')][0]
    list_offseqs  = [readline.strip('\n').replace(' ', '').split(',') for readline in open('%s/offseq.txt' % work_dir, 'r')][0]

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
                outfile.write('%s-off%s,%s,%s\n' % (pegID, (i+1), offseq, ','.join(dict_target[pegID])))
        outfile.close()
#def END: make_input


def off_target_test (work_dir):

    df              = pd.read_csv('%s/input_off.csv' % work_dir)
    df['Offtarget_Score'] = calculate_deepprime_score(sBASE_DIR, df, 'PE2_offtarget')
    sorted_df = df.sort_values(by=['Offtarget_Score'], ascending=False)
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


def mp_processor(list_sParameters):
    analysistag, work_dir, df_sSubSplits, dict_params = list_sParameters

    outdir   = '%s/output' % work_dir
    os.makedirs(outdir, exist_ok=True)

    list_output = []
    for idx in df_sSubSplits.index:
        df_temp         = df_sSubSplits.loc[idx]
        sID             = df_temp[0]
        Ref_seq         = df_temp[1].upper()
        ED_seq          = df_temp[2].upper()
        sAlt            = df_temp[3]
        result          = deep_prime(analysistag, sID, Ref_seq, ED_seq, sAlt, dict_params, outdir)
        list_output += result

        # print('Processing: %s' % sID)
        progress = round(100 * (idx / len(df_sSubSplits)), 1)
        print('Processing: %d%% - %s' % (progress, sID))
    # loop END: idx

    columns = ['ID', 'Total_pegRNAs', 'Average_DP_score', 'Average_Top4_DP_score',
               'Over30_pegRNAs', 'Over20_pegRNAs', 'Over10_pegRNAs', 'Over5_pegRNAs']
    header  = ['id', 'guide', 'pbsrtt', 'wtseq', 'edseq', 'spacer_top', 'spacer_bot', 'ext_top', 'ext_bot', 'dpcas', 'dpprime', 'pbslen', 'rttlen', 'edittype', 'zscore', 'percentile', 'plotpath']
    df_stat     = pd.DataFrame(list_output, columns=header)
    sorted_stat = df_stat.sort_values(by=['dpprime'], ascending=False)
    sorted_stat.to_csv('%s/Results_forHTML.csv' % work_dir, index=False, sep=',')

    print('Done')
# def END: mp_processor


def deep_prime(analysistag, sID, Ref_seq, ED_seq, sAlt, dict_params, outdir):

    nAltIndex   = ALT_INDEX
    pbs_range   = dict_params['PBS_range']
    rtt_max     = dict_params['RTT_max']
    pe_system   = dict_params['PE_system']
    edit_type   = sAlt[:-1]
    edit_len    = int(sAlt[-1])
    top         = 5

    results    = '%s/results' % outdir
    plots      = '%s/matplot' % outdir
    os.makedirs(results, exist_ok=True)
    os.makedirs(plots, exist_ok=True)

    plotphp    = '/DeepPrime/data/%s/output/matplot' % analysistag #for php

    ## FeatureExtraction Class
    cFeat          = FeatureExtraction()
    cFeat.input_id = sID
    cFeat.get_input(Ref_seq, ED_seq, edit_type, edit_len)

    cFeat.get_sAltNotation(nAltIndex)
    cFeat.get_all_RT_PBS(nAltIndex, nMinPBS=pbs_range[0] - 1, nMaxPBS=pbs_range[1], nMaxRT=rtt_max, pe_system=pe_system)
    cFeat.make_rt_pbs_combinations()
    cFeat.determine_seqs()
    cFeat.determine_secondary_structure()

    df = cFeat.make_output_df(bTEST)

    if len(df) == 0:  # Empty DataFrame = No PAM found
        return [cFeat.input_id, 0, 0, 0, 0, 0, 0, 0]

    else:
        list_Guide30 = [WT74[:30] for WT74 in df['WT74_On']]

        df['DeepSpCas9_score'] = calculate_DeepSpCas9_score(sBASE_DIR, list_Guide30)
        df['DeepPrime_score'], df['Zscores']  = calculate_deepprime_score(sBASE_DIR, df, pe_system)

        sorted_df       = df.sort_values(by=['DeepPrime_score'], ascending=False)

        df_top_pegs     = sorted_df[:top]
        dict_top        = df_top_pegs.to_dict('index')

        list_scores     = sorted(df['Zscores'])
        list_minmax     = minmax_scale(list_scores)

        list_topdesign  = get_oligos(plotphp, list_scores, dict_top)

        header          = 'id,guide,pbsrtt,wtseq,edseq,spacer_top,spacer_bot,ext_top,ext_bot,dpcas,dpprime,pbslen,rttlen,edittype,zscore,percentile,plotpath'
        outf = open('%s/%s.top_designs.csv' % (results, sID), 'w')
        outf.write('%s\n' % header)
        for data in list_topdesign:
            out = ','.join([str(col) for col in data])
            outf.write('%s\n' % out)
        #loop END: id
        outf.close()

        rank_plots (plots, list_scores, list_minmax, list_topdesign)

        ## Save result file
        # df.to_feather('%s/%s.feather' % (sOut_DIR, sID))
        sorted_df.to_csv('%s/%s.csv' % (results, sID), index=False)
        '''
        dp_score    = df.DeepPrime_score
        tot_pegRNAs = len(dp_score)
        ave_DPScore = np.mean(dp_score)
        ave_DPStop4 = np.mean(dp_score.sort_values(ascending=False).head(4))

        DPStop4     = df.sort_values(by='DeepPrime_score', ascending=False).head(4)

        over_5      = dp_score[dp_score >= 5]
        over_10     = over_5[over_5 >= 10]
        over_20     = over_10[over_10 >= 20]
        over_30     = over_20[over_20 >= 30]

        list_stat   = [cFeat.input_id, tot_pegRNAs, ave_DPScore, ave_DPStop4, len(over_30), len(over_20), len(over_10),
                     len(over_5)]
        '''
        return list_topdesign
    # if END:


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


def svg_clinvar (clinvarplot, pngfile, wtseq, edseq, editinfo):
    # Set window size
    svg_width    = 1500
    svg_height   = 300

    fontsize     = 20
    text_height  = 0
    text_height2 = 30

    list_script  = []
    list_script.append('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%d" height="%d">\n'
                        % (svg_width, svg_height))
    list_script.append('<rect width="%.2f" height="%.2f" style="fill:#ffffff;"/>\n'
                        % (svg_width, svg_height))
    list_script.append('<g transform="translate(%.3f,%.3f)" font-family="Arial">\n'
                        % (20, 40))

    #list_script.append('<text x="%.2f" y="%.2f" font-size="%.2f" text-anchor="middle" style="fill:%s;">%s</text>\n'
    #                        % (nText_X, nText_Y, nFontSize, sTextColor, sSegmentSeq[i]))

    list_script.append('<text x="10" y="%s" font-size="%s" style="font: font-weight:bold;">WTSeq: %s</text>\n'
                        % (text_height, fontsize, wtseq))
    list_script.append('<text x="10" y="%s" font-size="%s" style="font-weight:bold;">EDSeq: %s</text>\n'
                        % (text_height2, fontsize, edseq))

    list_script.append('</g>\n')
    list_script.append('</svg>\n')

    out      = ''.join(list_script)
    outf     = open(clinvarplot, 'w')
    outf.write(out)
    outf.close()

    ## [3] convert .svg into .png file
    svg2png(clinvarplot, pngfile)
#def END: plot_SS_on_SVG_compiled





def svg2png(infile, outpng):
    list_script = []
    list_script.append('java -jar ')
    list_script.append('/home/hkim/bin/batik-1.16/batik-rasterizer-1.16.jar ')
    list_script.append('%s -d %s > /dev/null' % (infile, outpng))

    script = ''.join(list_script)

    os.system(script)
#def END: svg2png


def main():
    print('main')
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
