from abf_reader import abf_reader as abr
from abf_reader import abf_epsd_cnstrctr as epsd_cnst
from abf_reader import dac_waveform
import numpy as np
from numpy import r_
import pdb
from pylab import plot, figure, plt
from glob import glob

def bwfiltfilt(data, sample_rate, cuttoff_freq):
    from scipy import signal
    import math
    # design filter
    norm_pass = 2*math.pi*cuttoff_freq/sample_rate
    norm_stop = 2*norm_pass
    (N, Wn) = signal.buttord(wp=norm_pass, ws=norm_stop, gpass=2, gstop=30, analog=0)
    (b, a) = signal.butter(N, Wn, btype='low', analog=0, output='ba')
    filtd = signal.filtfilt(b, a, data)
    return filtd

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


def main(abf_glob, chan_num, dac_num, epch_num, lpad_msec, rpad_msec, **kwds):
    from numpy import zeros, transpose
    abr_list = [abr(abf) for abf in abf_glob]

    # construct figure:
    from matplotlib.gridspec import GridSpec as gs
    num_plt_cols = len(abr_list)
    num_plt_rows = 1
    fig = plt.figure()
    mygs = gs(num_plt_rows,
              num_plt_cols,
              bottom = 0.3,
              top = 0.6)
    if 'ylim' in kwds.keys():
        ylim = kwds.pop('ylim')
    else:
        ylim = (-10,5)

    # loop over files in glob
    for col_num, abf in enumerate(abr_list):
        
        # get dac
        dd = dac_waveform(abf, dac_num)

        # loop over episodes in dac, filter each episode
        # first pre-index array
        num_epsds = abf._num_episodes
        epch = dd.get_epoch(epch_num)
        
        # this is dumb, have to deal with starts and ends
        sample_rate = abf.sample_rate()
        lpad_dp = sample_rate * lpad_msec / 1000
        rpad_dp  = sample_rate * rpad_msec / 1000
        padded_swp_len = epch._len() + lpad_dp + rpad_dp
        data = zeros((padded_swp_len, num_epsds))

        start = epch._dp_start() - epch.get_dp_pad() - lpad_dp
        end = epch._dp_end() - epch.get_dp_pad() + rpad_dp

        cntr = 0
        while 1:
            try:
                epsd_data = dd()
                filtd = bwfiltfilt(epsd_data[:,chan_num], sample_rate, 100)
                data[:,cntr] = filtd[start:end]
                dd.next_episode()
                cntr+=1
            except StopIteration:
                break

        # make axis
        ax = fig.add_subplot(mygs[0,col_num])
        
        # plot all episodes on one axis
        for epsd_num, epsd_data in enumerate(transpose(data)):
            ax.plot(epsd_data, **kwds)
        
        # set ylims
        ax.set_ylim(ylim)
    plt.show()

# debuggy stuff
# import os
# root_lab_dir = os.environ.get("LABDIR")
# from file_selector import same_size
# exp_path = 'B67_peptide/2012_03_08/binary_data/'
# fname = '2012_03_08_0007.abf'
# smpl_path = os.path.join(root_lab_dir, exp_path, fname)
# abf_file_list = same_size(os.path.dirname(smpl_path), fname)

if __name__=='__main__':
    # select first as representative abf
    exmplr = abr(abf_file_list[0])

    # select chan name
    exmplr.chan_names()
    chan_num = input('which channel\n')

    # select dac num
    dac_num = input('which dac_num\n')

    # select epoch num
    epoch_num = input('which epoch num\n')

    # get lpad msec
    lpad_msec = input('left pad milli seconds\n')

    # get rpad msec
    rpad_msec = input('right pad milli seconds\n')

    # chan_num = 1
    # dac_num = 0
    # epoch_num = 2
    # lpad_msec = 50
    # rpad_msec = 100

    main(abf_file_list,
         chan_num,
         dac_num,
         epoch_num,
         lpad_msec,
         rpad_msec,
         color = 'blue', ylim = (-20,20))
