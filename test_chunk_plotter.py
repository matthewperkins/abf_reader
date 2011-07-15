from abf_reader import abf_reader as abr
from pylab import *
import pdb
import chunker
import os

tabr = abr(os.path.join( os.environ['HOME'], 'Documents',
                         'weiss_lab', 'B8_experiments', '2011_06_30',
                         'binary_data', '2011_06_30_0011.abf'))

test_chunker = chunker.abf_chunker(tabr)
test_chunker.set_range((0,160*5000))
test_chunker.make_chunk_size(2**22)
phi = 1.61803399

class abf_chunker_plotter(object):
    def __init__(self, _abf_chunker):
        self.num_chunks = _abf_chunker.num_chunk_rows
        self._abf_chunker = _abf_chunker
        self.frct_chunks = self._abf_chunker.frct_chunks
        self.set_fig_size_inches((8,3))
        self._ylim = (-85,20)

    def set_ylim(self, rng):
        self._ylim = rng

    def set_fig_size_inches(self, size):
        self.indiv_fig_width = float(size[0]) / (self.frct_chunks)
        self.indiv_fig_height = float(size[1])

    def plot(self, chan, dpi):
        for i, d in enumerate(self._abf_chunker):
            data = d[:,chan]
            fig = plt.figure()
            fig.patch.set_alpha(0)
            ax = plt.axes([0,0,1,1])
            ax.patch.set_alpha(0)
            ax.plot(data, linewidth = 0.2, color = '#0F0B3B')
            ax.set_ylim(self._ylim)
            ax.set_xlim((0,len(d)))
            ax.set_axis_off()
            fig.add_axes(ax)
            # want to 
            width_frct = float(len(d))/self._abf_chunker.chunk_row_size
            fig.set_size_inches(( self.indiv_fig_width * width_frct, self.indiv_fig_height))
            filnum = "%03d" % (i)
            fig.savefig('tmp' + filnum + '.png', dpi = dpi)
            plt.close()

chnk_plter = abf_chunker_plotter(test_chunker)
chnk_plter.set_ylim((-70, 20))
chnk_plter.set_fig_size_inches((8,2))
tabr.chan_names()
chan_no = raw_input('give me the chan number you want')
chnk_plter.plot(int(chan_no),300)


# then use imagemagick to stich:
# montage -tile x1 -background transparent -mode Concatenate `ls *.png` out.png

