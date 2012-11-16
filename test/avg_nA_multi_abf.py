from abf_reader import abf_reader as abr
from abf_reader import abf_epsd_cnstrctr as epsd_cnst
from abf_reader import dac_waveform
import numpy as np
from numpy import r_
import pdb
from pylab import plot, figure, plt
from glob import glob

def avg_nA(abf_list, dac_num, epch_num, strt_dp, end_dp):
    dac_num = dac_num
    epch_strt_dp = strt_dp
    epch_end_dp = end_dp

    # pre index two array, one for voltage, other for current.
    # arrays will be [file_num, epsd_num] in order

    # rows is num files
    nrows = len(abf_list)

    # cols is num epchs, assume first file similar to list
    smpl_abf = abf_list[0]
    ncols = smpl_abf._num_episodes

    # make array
    nA = np.zeros((nrows, ncols), dtype = np.float32)
    Vm = np.zeros((nrows, ncols), dtype = np.float32)

    for fnum, abr_f in enumerate(abf_list):
        dd = dac_waveform(abr_f, dac_num)
        epch = dd.get_epoch(epch_num)
        for epsd_num, epch_epsd in enumerate(epch):
            vmem = np.average(epch_epsd()[epch_strt_dp:epch_end_dp,0])
            imem = np.average(epch_epsd()[epch_strt_dp:epch_end_dp,1])
            nA[fnum, epsd_num] = imem
            Vm[fnum, epsd_num] = vmem
    return (nA, Vm)

def main(abf_glob, start_dp=2500, end_dp=3000):
    abr_list = [abr(abf) for abf in abf_glob]
    nA, Vm = avg_nA(abr_list, 0, 2, start_dp, end_dp)
    for col in range(nA.shape[1]):
        plot(nA[:,col], 'o')
    fig = plt.gcf()
    fig_name = "%s_%s.%s" % (start_dp, end_dp, 'pdf')
    fig.savefig(fig_name)
    plt.show()

# debuggy stuff
import os
root_lab_dir = os.getenv("LABDIR")
from file_selector import same_size
exp_path = 'B67_peptide/2012_03_08/binary_data/'
fname = '2012_03_08_0007.abf'
smpl_path = os.path.join(root_lab_dir, exp_path, fname)
abf_file_list = same_size(os.path.dirname(smpl_path), fname)
if __name__ == '__main__':
    main(abf_file_list)

