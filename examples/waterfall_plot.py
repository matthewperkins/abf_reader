from abf_reader import *
import matplotlib.pyplot as plt
from abf_epoch import epoch, waterfall
import matplotlib.gridspec as gridspec
from scale_bars import *

if __name__=='__main__':
    import os
    labdir = os.environ.get("LABDIR")
    
    # change some legend plot stuff
    plt.rcParams.update(\
        {'legend.frameon':False,
         'legend.numpoints':1,
         'legend.loc':'upper left',
         'legend.fontsize':12})

    # organize info for files / conditions / plot labels into dict
    files_dict = [{'abf':abf_reader('2012_08_24_0000_cntrl.abf'),
             'color':'black',
             'label':'control'},
            {'abf':abf_reader('2012_08_24_0004_frf.abf'),
             'color':'red',
             'label':'frf'}]

    # pre plot grid
    plt.figure()
    gs = gridspec.GridSpec(2, 1, hspace = 0, wspace = 0,
                           left = 0, right = 1.)

    # plot in loop
    xlims = []
    for i, fdict in enumerate(files_dict):
        plt.subplot(gs[i,0])
        epch = epoch(fdict.pop('abf'), 1, 1)
        epch.set_pading(left = 1000, right = 2000)
        wf = waterfall(epch)
        wf.set_range() # this will change the xlim
        xlims.append(wf._xlim)
        for i, (xs,ys) in enumerate(wf):
            plt.plot(xs,ys, linewidth = 0.5, **fdict)
            if i==0:
                fdict.pop('label')
        plt.show()

    # post plotting adjust axis so grids are comparable, maybe can do
    # this with some shared axis object? Also hide axis ticks, labels
    # and spines.
    xmin = min([x[0] for x in xlims])
    xmax = max([x[1] for x in xlims])
    for ax in plt.gcf().axes:
        ax.set_xlim((xmin, xmax)) 
        ax.set_ylim((-60,30))
        ax.legend()

        # add a scale bar
        add_scalebar(ax, matchx = False, matchy = False,
                     sizex = 1, sizey = 30,
                     labelx = '1 sec', labely = '30 mV',
                     bbox_to_anchor = (0.3,0.4),
                     sep = 1, pad = 1,
                     bbox_transform = ax.transAxes, borderpad=0)
        plt.draw()

