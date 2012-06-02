from abf_reader import abf_reader as abr
from abf_reader import abf_epsd_cnstrctr as epsd_cnst
from abf_reader import dac_waveform
import numpy as np
from numpy import r_
import pdb
from pylab import plot, figure, plt

import os
lab_dir = os.environ.get("LABDIR")
t = abr(os.path.join(lab_dir, 'B44_B48_experiments' ,
                     '2012_03_02' ,
                     'binary_data', '2012_03_02_0008.abf'))

dac_num = 0
c = epsd_cnst(t)
elist = []
while 1:
    try:
        elist.append(c.next())
    except StopIteration:
        break

for epsd in elist:
    plot(epsd()[:,1], color = 'black', linewidth = 0.5)

figure()
dd = dac_waveform(t, 0)
while 1:
    try:
        dd.next_episode()
        plot(dd(cmd = 'level'))
    except StopIteration:
        dd.rewind_episode()
        break

epch_strt_dp = 300
epch_end_dp = 1000

epch1 = dd.waveforms[1]
figure()
for epsd_epch in epch1:
        vm = np.average(epsd_epch()[epch_strt_dp:epch_end_dp,0])
        im = np.average(epsd_epch()[epch_strt_dp:epch_end_dp,1])
        plot(epsd_epch()[epch_strt_dp:epch_end_dp,0], color = 'black')
        plot(epsd_epch()[epch_strt_dp:epch_end_dp,1], color = 'red')
plt.gca().set_ylim((-100,40))
plt.show()      

### finally FIGURED OUT HOW THE PRE AND POST HOLDS ARE DETERMINED,
### THEY ARE 1/64 OF THE EPISODE LENGTH
