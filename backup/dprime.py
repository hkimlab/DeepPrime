# Data for normalization should be in models directory

import os, sys
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import torch.nn as nn
from glob import glob
from src.utils import seq_concat, select_cols



class GeneInteractionModel(nn.Module):

    def __init__(self, hidden_size, num_layers, num_features=24, dropout=0.1):
        super(GeneInteractionModel, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers

        self.c1 = nn.Sequential(
            nn.Conv2d(in_channels=4, out_channels=128, kernel_size=(2, 3), stride=1, padding=(0, 1)),
            nn.BatchNorm2d(128),
            nn.GELU(),
        )
        self.c2 = nn.Sequential(
            nn.Conv1d(in_channels=128, out_channels=108, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(108),
            nn.GELU(),
            nn.AvgPool1d(kernel_size=2, stride=2),

            nn.Conv1d(in_channels=108, out_channels=108, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(108),
            nn.GELU(),
            nn.AvgPool1d(kernel_size=2, stride=2),

            nn.Conv1d(in_channels=108, out_channels=128, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(128),
            nn.GELU(),
            nn.AvgPool1d(kernel_size=2, stride=2),
        )

        self.r = nn.GRU(128, hidden_size, num_layers, batch_first=True, bidirectional=True)

        self.s = nn.Linear(2 * hidden_size, 12, bias=False)

        self.d = nn.Sequential(
            nn.Linear(num_features, 96, bias=False),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(96, 64, bias=False),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 128, bias=False)
        )

        self.head = nn.Sequential(
            nn.BatchNorm1d(140),
            nn.Dropout(dropout),
            nn.Linear(140, 1, bias=True),
        )

    def forward(self, g, x):
        g = torch.squeeze(self.c1(g), 2)
        g = self.c2(g)
        g, _ = self.r(torch.transpose(g, 1, 2))
        g = self.s(g[:, -1, :])

        x = self.d(x)

        out = self.head(torch.cat((g, x), dim=1))

        return F.softplus(out)


def calculate_deepprime_score(sBase_DIR, df_input, pe_system='PE2'):

    os.environ['CUDA_VISIBLE_DEVICES']='0'
    device     = 'cuda' if torch.cuda.is_available() else 'cpu'

    model_dir  = '%s/bin/models' % sBase_DIR
    dict_model = { 'PE2':                   'DeepPrime_base',
                   'PE2_offtarget':         'DeepPrime_off',

                   'PE2_Conv':              'DeepPrime_var/DP_variant_293T_PE2_Conv_220428',
                   'NRCH-PE2_HEK293T':      'DeepPrime_var/DP_variant_293T_NRCH_PE2_Opti_220428',

                   'PE2_Opti_HCT116':       'DeepPrime_var/DP_variant_HCT116_PE2_Opti_220428',
                   'PE2_Opti_MDA':          'DeepPrime_var/DP_variant_MDA_PE2_Opti_220428',

                   'PE2max_Opti_HEK239T':   'DeepPrime_var/DP_variant_293T_PE2max_Opti_220428',
                   'PE2max_Opti_Hela':      'DeepPrime_var/DP_variant_HeLa_PE2max_Opti_220815',
                   'PE2max_Opti_A549':      'DeepPrime_var/DP_variant_A549_PE2max_Opti_221114',
                   'PE2max_Opti_DLD1':      'DeepPrime_var/DP_variant_DLD1_PE2max_Opti_221114',

                   'NRCH-PE2max_Opti_HEK293T':      'DeepPrime_var/DP_variant_293T_NRCH-PE2max_Opti_220815',

                   'PE4max_Opti_HEK293T':           'DeepPrime_var/DP_variant_293T_PE4max_Opti_220728',
                   'PE4max_Opti_A549':              'DeepPrime_var/DP_variant_A549_PE4max_Opti_220728',
                   'PE4max_Opti_DLD1':              'DeepPrime_var/DP_variant_DLD1_PE4max_Opti_220728',

                   'NRCH-PE4max_Opti_DLD1':         'DeepPrime_var/DP_variant_DLD1_NRCHPE4max_Opti_220728',
                   'NRCH-PE4max_Opti_NIH':          'DeepPrime_var/DP_variant_NIH_NRCHPE4max_Opti_220815',

                   'PE2max_epegRNA_Opti_HEK293T':   'DeepPrime_var/DP_variant_293T_PE2max_epegRNA_Opti_220428',
                   'PE2max_epegRNA_Opti_A549':      'DeepPrime_var/DP_variant_A549_PE2max_epegRNA_Opti_220428',

                   'PE4max_epegRNA_Opti_HEK293T':   'DeepPrime_var/DP_variant_293T_PE4max_epegRNA_Opti_220428',
                   'PE4max_epegRNA_Opti_A549':      'DeepPrime_var/DP_variant_A549_PE4max_epegRNA_Opti_220428',
                  }


    model_type  = dict_model[pe_system]
    mean        = pd.read_csv('%s/%s/mean.csv' % (model_dir, model_type), header=None, index_col=0, squeeze=True) # Train set mean (made with preprocessing.py)
    std         = pd.read_csv('%s/%s/std.csv'  % (model_dir, model_type), header=None, index_col=0, squeeze=True) # Train set std (made with preprocessing.py)

    dpmean = float(open('%s/%s/dp_mean.csv' % (model_dir, model_type), 'r').readline())  #  DP scores from test set
    dpstd  = float(open('%s/%s/dp_std.csv' % (model_dir, model_type), 'r').readline())  #  DP scores from test set

    test_feats, _  = select_cols(df_input)

    if model_type == 'DeepPrime_off':
        g_test = seq_concat(df_input, col1='WT74_ref', col2='Edited74_On')
    else:
        g_test = seq_concat(df_input)


    x_test      = (test_feats - mean) / std

    g_test      = torch.tensor(g_test, dtype=torch.float32, device=device)
    x_test      = torch.tensor(x_test.to_numpy(), dtype=torch.float32, device=device)

    models      = [m_files for m_files in glob('%s/%s/*.pt' % (model_dir, model_type))]
    preds       = []

    for m in models:
        model = GeneInteractionModel(hidden_size=128, num_layers=1).to(device)
        model.load_state_dict(torch.load(m))
        model.eval()
        with torch.no_grad():
            g, x    = g_test, x_test
            g       = g.permute((0, 3, 1, 2))
            pred    = model(g, x).detach().cpu().numpy()
        preds.append(pred)
    #loop END: m

    # AVERAGE PREDICTIONS
    preds = np.squeeze(np.array(preds))
    preds = np.mean(preds, axis=0)
    preds = np.exp(preds) - 1

    zscores = []
    for pred in preds:
        zscores.append((pred - dpmean) / dpstd)

    return preds,zscores

# # SAVE RESULTS
# preds = pd.DataFrame(preds, columns=['Predicted_PE_efficiency'])
# preds.to_csv('prediction.csv', index=False)

