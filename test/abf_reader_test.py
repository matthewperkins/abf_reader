from abf_reader import abf_reader, dac_waveform
from pylab import *

import os
import sys
import pdb

abf_path = os.path.join(os.getenv("LABDIR"),
                         'matts_axon18',
                         'matts_axon',
                         'test')
abf = os.path.join(abf_path, 'test.abf')
abr = abf_reader(abf)

# make a waveform object from the epochs defined in dac 0 
dac_wvf = dac_waveform(abr, 0)
epch1 = dac_wvf.get_epoch(2)

# for selecting channels
chn_no = 0
# or chn_slc = slice(strt_chan, stop_chan)

sweep_slc = slice(None, None)
for sweep in epch1:
    # get data with a call, select chan and 
    d = sweep()[sweep_slc, chn_no]

    #plot
    plot(d)

# another way to do this is by
while 1:
    try:
        dac_wvf.next_episode()
        d = dac_wvf.get_epoch(2)()[sweep_slc, chn_no]
        plot(d)
    except StopIteration:
        break

### plot all sweeps
dac_wvf.rewind_episode()    
sweep_slc = slice(-40000, -20000)
plt.figure()
while 1:
    try:
        dac_wvf.next_episode()
        pdb.set_trace()
        d = dac_wvf()[sweep_slc, chn_no]
        plot(d)
    except StopIteration:
        break
    
plt.show()
    

