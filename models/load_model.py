import os, sys, glob
import inspect

from models import DeepSpCas9, DeepPrime

'''
Model finder repo
Get your models from here!
'''

def load_deepspcas9():
    '''
    Load DeepSpCas9 model
    

    '''

    print('get model: DeepSpCas9')

    model_dir = inspect.getfile(DeepSpCas9).replace('/__init__.py', '').replace('\\__init__.py', '')

    return model_dir


def load_deepprime(model_id='PE2', cell_type='HEK293T'):
    '''
    model_id: PE2, PE2max, PE4max, PE2max-e, PE4max-e, NRCH_PE2, NRCH_PE2max, NRCH_PE4max
    cell_rtpe: HEK293T, A549, DLD1, HCT116, HeLa, MDA-MB-231, NIH3T3

    >>> from genet_models import load_model
    >>> model_dir, model_type = load_model('PE2max', 'HEK293T')

    '''

    print('get model: %s - %s' % (model_id, cell_type))

    model_dir = inspect.getfile(DeepPrime).replace('/__init__.py', '').replace('\\__init__.py', '')

    dict_models = {
        
        'HEK293T': {
            'PE2'        : 'DeepPrime_base',
            'NRCH_PE2'   : 'DP_variant_293T_NRCH_PE2_Opti_220428',
            'NRCH_PE2max': 'DP_variant_293T_NRCH-PE2max_Opti_220815',
            'PE2max'     : 'DP_variant_293T_PE2max_Opti_220428',
            'PE2max-e'   : 'DP_variant_293T_PE2max_epegRNA_Opti_220428',
            'PE4max'     : 'DP_variant_293T_PE4max_Opti_220728',
            'PE4max-e'   : 'DP_variant_293T_PE4max_epegRNA_Opti_220428',
        },

        'A549': {
            'PE2max'     : 'DP_variant_A549_PE2max_Opti_221114',
            'PE2max-e'   : 'DP_variant_A549_PE2max_epegRNA_Opti_220428',
            'PE4max'     : 'DP_variant_A549_PE4max_Opti_220728',
            'PE4max-e'   : 'DP_variant_A549_PE4max_epegRNA_Opti_220428',

        },
        
        'DLD1': {
            'NRCH_PE4max': 'DP_variant_DLD1_NRCHPE4max_Opti_220728',
            'PE2max'     : 'DP_variant_DLD1_PE2max_Opti_221114',
            'PE4max'     : 'DP_variant_DLD1_PE4max_Opti_220728',
        },

        'HCT116': {
            'PE2'        : 'DP_variant_HCT116_PE2_Opti_220428',
        },
        
        'HeLa': {
            'PE2max'     : 'DP_variant_HeLa_PE2max_Opti_220815',
        },
        
        'MDA-MB-231': {
            'PE2'        : 'DP_variant_MDA_PE2_Opti_220428',
        },
        
        'NIH3T3': {
            'NRCH_PE4max': 'DP_variant_NIH_NRCHPE4max_Opti_220815',
        },

    }


    try:
        model_type = dict_models[cell_type][model_id]
    except:
        print('Not available Prime Editor')
        sys.exit()

    
    return model_dir, model_type