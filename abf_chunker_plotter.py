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
    elif 'cell_anchors' and 'cell_spans' in kwargs.keys():
        anchrs = kwargs.pop('cell_anchors')
        spns = kwargs.pop('cell_spans')
        cell_lms = [(anchr, anchr+spn) for anchr, spn in zip(anchrs, spns)]
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

    # lets set some defaults:
    options={'cells':range(abrf.num_chans()),
             'neurgrms':[],
             'dpi':300,
             'srange':(0.,100.),
             'lnwdth':0.3,
             'cell_lms':[(-80,30)]*abrf.num_chans(),
             'neurgrm_lms':[(-3,3)*0],
             'fig_sz':(10,8),
             'filetype':'png'}

    for key in kwargs.keys():
        if key=='cells':
            if 'cell_lms' not in kwargs.keys():
                options['cell_lms'] = [(-80,30)]*len(kwargs['cells'])
        if key=='neurgrms':
            if 'neurgrm_lms' not in kwargs.keys():
                options['neurgrm_lms'] = [(-3,3)]*len(kwargs['cells'])
        options[key] = kwargs[key]
    # give chunker an abr 
    abf_chunker = chunker.abf_chunker(abrf)

    # set up chunker plotter
    chnk_plter = abf_chunker_plotter(abf_chunker)

    # set line width
    slinwidth(chnk_plter, **options)
    
    # select cells
    if 'cells' in options.keys():
        cells = options.pop('cells')
    else:
        cells = scells(abrf)

    # select neurgrams
    if 'neurgrms' in options.keys():
        neurgrms = options.pop('neurgrms')
    else:
        neurgrms = sneurgrms(abrf)

    # set yrng
    syrng(chnk_plter, cells, neurgrms, **options)

    # set xrange in seconds
    sxrange(abf_chunker, **options)

    # set fig size
    sfigsz(chnk_plter, **options)

    # get dpi
    if 'dpi' in options.keys():
        dpi = options.pop('dpi')
    else:
        dpi = input('what dpi?\n')

    chnk_plter.plot(dpi, cells, neurgrms, **options)
    # changes:
    
if __name__ == '__main__':
    main(sys.argv[1])

# debug stuff
    
