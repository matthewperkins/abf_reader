#!/usr/bin/env python
import pdb
from abf_reader import abf_reader as abr
from pylab import *
from chunk_plotter import abf_chunker_plotter
import os
import chunker
import sys

def sxrange(abf_chunker, **kwargs):
    # set range in seconds
    if 'srange' in kwargs.keys():
        srange = kwargs.pop('srange')
    else:
        srange = input('give me a range tuple of floats, for the time range in seconds\n')
    abf_chunker.set_range(srange, seconds = True)

def slinwidth(chnk_plter, **kwargs):
    if 'lnwdth' in kwargs.keys():
        lnwdth = kwargs.pop('lnwdth')
    else:
        # set line width
        lnwdth = input('give me a linewidth, 0.4 works okay\n')
    chnk_plter.set_linewidth(lnwdth)

def scells(abrf):
    # list channel names
    abrf.chan_names()
    # get cells to plot from promt
    cells = input('give me a list of intra cells\n')
    return cells

def sneurgrms(abrf):
    # list channel names
    abrf.chan_names()
    # get neurgrams to plot from prompt
    neurgrms = input('give me a list of Neurgrams\n')
    return neurgrms

def syrng(chnk_plter, cells, neurgrms, **kwargs):
    if 'cell_lms' in kwargs.keys():
        cell_lms = kwargs.pop('cell_lms')
    else:
        special_cell = input('special cells yrange? 1 for yes\n')
        if special_cell==1:
            cell_lms = input('give a list of tuples, one for each cell\n')
        else:
            cell_lms = [(-80,30)]*len(cells)
    if 'neurgrm_lms' in kwargs.keys():
        neurgrm_lms = kwargs.pop('neurgrm_lms')
    else:
        special_neurgrm = input('special neurgrams yrange? 1 for yes\n')
        if special_neurgrm==1:
            neurgrm_lms = input('give a list of tuples, one for each neurogram\n')
        else:
            neurgrm_lms = [(-3,3)]*len(neurgrms)
    chnk_plter.set_ylim(neurgrms = neurgrm_lms, cells = cell_lms)

def sfigsz(chnk_plter, **kwargs):
    if 'fig_sz' in kwargs.keys():
        fig_sz = kwargs.pop('fig_sz')
    else:
        fig_sz = input('give me a tuple for the figure size l x w (inches)\n')

    chnk_plter.set_fig_size_inches(fig_sz)

def main(abf_path, **kwargs):
    # make abf reader instance
    abrf = abr(abf_path)
    # print channel names, for sanity
    abrf.chan_names()
    print('\n')
    # set up the abf_chunker
    # give chunker an abr 
    abf_chunker = chunker.abf_chunker(abrf)

    # set up chunker plotter
    chnk_plter = abf_chunker_plotter(abf_chunker)

    # set line width
    slinwidth(chnk_plter, **kwargs)
    
    # select cells
    if 'cells' in kwargs.keys():
        cells = kwargs.pop('cells')
    else:
        cells = scell(abrf)

    # select neurgrams
    if 'neurgrms' in kwargs.keys():
        neurgrms = kwargs.pop('neurgrms')
    else:
        neurgrms = sneurgrms(abrf)

    # set yrng
    syrng(chnk_plter, cells, neurgrms, **kwargs)

    # set xrange in seconds
    sxrange(abf_chunker, **kwargs)

    # set fig size
    sfigsz(chnk_plter, **kwargs)

    # get dpi
    if 'dpi' in kwargs.keys():
        dpi = kwargs.pop('dpi')
    else:
        dpi = input('what dpi?\n')

    chnk_plter.plot(dpi, cells, neurgrms, **kwargs)
    # changes:
    
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

# debug stuff
    