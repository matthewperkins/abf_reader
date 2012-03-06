from abf_reader import abf_reader as abr
from abf_reader import abf_epsd_cnstrctr as epsd_cnst
from abf_reader import dac_waveform
import numpy as np
from file_selector import root_lab_dir as rld
from numpy import r_
import pdb
from pylab import plot, figure, plt

import os
t = abr(os.path.join(os.path.expanduser('~'), 'Documents/weiss_lab/B44_B48_experiments/2012_03_02/binary_data/2012_03_02_0008.abf'))

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

dd = dac_waveform(t, 0)
        
epch1 = dd.waveforms[1]
figure()
while 1:
    try:
        epch1._next_episode()
        plot(epch1()[:,1], color = 'black')
        plot(epch1.level(), color = 'red')
    except StopIteration:
        break
plt.gca().set_ylim((-40,50))
plt.show()      

### finally FIGURED OUT HOW THE PRE AND POST HOLDS ARE DETERMINED,
### THEY ARE 1/64 OF THE EPISODE LENGTH
