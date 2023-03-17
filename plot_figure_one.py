"""
NOTE: units for "e_Mass" in "t1_mrf.txt" should be changed to "solMass"
      for compatibility with astropy.Table.read "ascii.cds" format

"""

# ===================== Import Packages ====================

import sys, os, pdb, glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import Table, join, MaskedColumn


# ===================== Define Functions ===================

def spt_coding(spt):

    """
    PURPOSE:    Translate spectral type (e.g., M7) into numerical value
                that can be used for plotting a histogram

                Scale is 0 at M0, -1 at K7, -8 at K0, -18 at G0
                (K8 is counted as M0)

    INPUT:      Numpy array of spectral types (str, masked)
                (masked values are unknonw spectral types)

    OUTPUT:     Numpy array of numerical spectral types (obj, masked)
                (returns -99. if unknown spectral type)
    
    """

    spt_num = np.empty(len(spt), dtype=object)
    for i, val in enumerate(spt):

        if np.ma.is_masked(val):
            spt_num[i] = -99.

        else:

            if val == 'M1+M2':
                val = 'M1.5'
            if val[0] == 'M':
                spt_num[i] = float(val[1:])
            if val[0] == 'K':
                spt_num[i] = float(val[1:]) - 8.
            if val[0] == 'G':
                spt_num[i] = float(val[1:]) - 18.
            if val[0] == 'F':
                spt_num[i] = float(val[1:]) - 28.
            if val[0] == 'A':
                spt_num[i] = float(val[1:]) - 38.
 
    return spt_num


# ========================== Code ==========================

#### LOAD IN TABLES FROM PAPER SUPPLEMENTAL MATERIAL
Table_Sample = Table.read('input/t1_mrf.txt', format='ascii.cds')
Table_Obs = Table.read('input/t2_mrf.txt', format='ascii.cds')
Table_All = join(Table_Sample, Table_Obs, join_type='outer')

### FLAG UNDETECTED SOURCES
Table_All.add_column(MaskedColumn(name='Detected', data=Table_All['FCont']/Table_All['e_FCont'] > 3.0))

### GET SPECTRAL TYPE CODING FOR PLOTTING (ALL, NON-DETECTIONS, UNOBSERVED)
Spt_All = spt_coding(Table_All['SpType'])
Spt_Ndt = spt_coding(Table_All['SpType'][np.where(~Table_All['Detected'])])
Spt_Obs = spt_coding(Table_All['SpType'][~Table_All['FCont'].mask])

### SETUP PLOT
mpl.rc('xtick', labelsize=15)
mpl.rc('ytick', labelsize=15)
mpl.rc('xtick.major', size=3, pad=7, width=1.5)
mpl.rc('ytick.major', size=3, pad=7, width=1.5)
mpl.rc('axes', linewidth=1.5)
mpl.rc('lines', markersize=5)
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
ax.set_xlabel('Spectral Type', fontsize=17)
ax.set_ylabel("Number of Sources", fontsize=17)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
nxticks    = np.arange(-33, 9, 1)
ticklabels = np.array(['A5', '', '', '', '', 'F0', '', '', '', '', 'F5', '', '', '', '', 'G0', '', '', '', '', 'G5', '', '', '', '',
                        'K0', '', '', '', '', 'K5', '', '', 'M0', '', '', 'M3', '', '', '', 'M7', ''])
ax.set_xticks(nxticks)
ax.set_xticklabels(ticklabels)
ax.set_xlim(np.min(nxticks), np.max(nxticks))
ax.set_ylim(0, 25)

### PLOT DIFFERENT SAMPLES
bins_spt = np.arange(-33, 9, 2)
ax.hist(Spt_All, bins_spt, align='left', lw=2, histtype='step', color='black', label='All Sources')
ax.hist(Spt_Obs, bins_spt, align='left', facecolor='lightblue', hatch='//', lw=2, histtype='stepfilled', label='Observed', edgecolor='black')
ax.hist(Spt_Ndt, bins_spt, align='left', color='red', lw=2, histtype='step', label='Undetected')
ax.legend(loc='upper left', prop={'size': 16}, edgecolor='black')

### SAVE FIGURE
fig.savefig('output-test/figure_one.pdf', bbox_inches='tight', dpi=100)
plt.close('all')


