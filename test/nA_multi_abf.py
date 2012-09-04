from abf_reader import abf_reader as abr
from abf_reader import abf_epsd_cnstrctr as epsd_cnst
from abf_reader import dac_waveform
import numpy as np
from numpy import r_
import pdb
from pylab import plot, figure, plt
from glob import glob

def many_nA(abf_list, dac_num, epch_num, epsd_num):
    dac_num = dac_num

    # pre index two array, one for voltage, other for current.
    # arrays will be [file_num, epch_dp] in order

    # rows is num files
    nrows = len(abf_list)

    # cols are each dp in epch_epsd
    smpl_abf = abf_list[0]
    dd = dac_waveform(smpl_abf, dac_num)
    epch = dd.waveforms[epch_num]
    epch._set_epsd_num(epsd_num)
    ncols = epch._len()

    # make array
    nA = np.zeros((nrows, ncols), dtype = np.float32)

    for fnum, abr_f in enumerate(abf_list):
        dd = dac_waveform(abr_f, dac_num)
        epch = dd.waveforms[epch_num]
        epch._set_epsd_num(epsd_num)
        nA[fnum,:] = epch()[:,1]
    return nA

def ad_hoc_plot(fat_array, **kwds):
    from matplotlib.gridspec import GridSpec as gs
    num_plt_cols = fat_array.shape[0]
    num_plt_rows = 1
    fig = plt.figure()
    mygs = gs(num_plt_rows,num_plt_cols, bottom = 0.3, top = 0.6)
    for row_num, row in enumerate(fat_array):
        ax = fig.add_subplot(mygs[0,row_num])
        ax.plot(row)
        if 'ylim' in kwds.keys():
            ax.set_ylim(kwds['ylim'])
    mygs.tight_layout(fig, pad = 0.1, w_pad = -2.5)


def main(abf_glob):
    abr_list = [abr(abf) for abf in abf_glob]

    # do user inter to get these values
    # assume first in list is representative
    exmplr = abr_list[0]

    # select dac
    print('active dacs')
    print(exmplr.actv_dacs())
    dac_num = input('which dac num')
    dd = dac_waveform(exmplr, dac_num)

    # select epoch
    print('active epochs')
    print(dd._actv_epchs)
    epch_num = input('which epch num')

    # select episode
    print('episodes')
    print(range(exmplr._num_episodes))
    epsd_num = input('which episode')
    
    # get data
    nA = many_nA(abr_list, dac_num, epch_num, epsd_num)

    #plot_data
    ad_hoc_plot(nA, ylim = (0,30))
    plt.show()

# debuggy stuff
import os
root_lab_dir = os.environ.get("LABDIR")
from file_selector import same_size
exp_path = 'B67_peptide/2012_03_08/binary_data/'
fname = '2012_03_08_0007.abf'
smpl_path = os.path.join(root_lab_dir, exp_path, fname)
abf_file_list = same_size(os.path.dirname(smpl_path), fname)
abf_file_a = np.array(abf_file_list)
if __name__=='__main__':
    main(abf_file_a)
