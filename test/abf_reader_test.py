from abf_reader import abf_reader
from pylab import *

import os
import sys
import glob

abfs_path = os.path.expanduser('~') + '/Downloads/testabfs'
abfs = glob.glob(abfs_path + os.sep + '*.abf')
abfs.sort()
anchr = abf_reader(abfs[0])
anchr_time = anchr.start_time()
for abf in abfs:
    abr = abf_reader(abf)
    abr.episodes_from_header()
    abr.epochs_from_header()

    #want to plot the data from the 2nd epond, (0 indexing)
    dpstart = abr.epochs[1]._dp_start
    dpend =  abr.epochs[1]._dp_start + abr.epochs[1]._dur_init
    abr.episode_data()

    #for each episode, take the average of the last 1000 data points at the
    #end of epoch 1(zero indexing), and average along the axis of
    #expisode(mean(1) give average along this axis)
    ys = abr.data[0:-1,dpend-1000:dpend]['adc_3'].mean(1)

    #want to plot these data according to the relative file start time, from
    #the anchr_time of the first file in the series.
    time_diff = abr.start_time() - anchr_time
    time_diff_seconds = time_diff.total_seconds()
    x = np.array(time_diff_seconds)
    #need to make to times the same shape as the y values:
    xs = np.tile(x, ys.shape)
    plot(xs,ys, 'ko')

