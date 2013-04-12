from abf_reader import abf_reader as abr
from pylab import *
from chunk_plotter import abf_chunker_plotter
import os
import chunker
import sys

def main(abf_path):

    # make abf reader instance
    abrf = abr(abf_path)

    # set up the abf_chunker
    # give chunker an abr 
    test_chunker = chunker.abf_chunker(abrf)

    # set range in seconds
    srange = input('give me a range tuple of floats')
    test_chunker.set_range(srange, seconds = True)

    # set up chunker plotter
    chnk_plter = abf_chunker_plotter(test_chunker)


    # list channel names
    abrf.chan_names()

    # get cells to plot from promt
    cells = input('give me a list of intra cells\n')
    cell_lims = [(-70,30)]*len(cells)

    # get neurgrams to plot from prompt
    neurgrms = input('give me a list of Neurgrams\n')
    neurgrm_lms = [(-300,300)]*len(neurgrms)

    # set ylims
    chnk_plter.set_ylim(neurgrms = neurgrm_lms, cells = cell_lims)

    # set fig size
    chnk_plter.set_fig_size_inches((10,8))
    chnk_plter.plot(96,cells, neurgrms)

if __name__ == '__main__':
    main(sys.argv[1])

