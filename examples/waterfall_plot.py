from abf_reader import *
import matplotlib.pyplot as plt
from abf_epoch import epoch, waterfall
import matplotlib.gridspec as gridspec
from scale_bars import *

default_sb_prop = {'matchx':False, 'matchy':False,
                   'sizex':1, 
                   'labelx':'1 sec',
                   'sizey':30,
                   'labely':'30 mV',
                   'bbox_to_anchor':(0.2,0.2),
                   'sep':1, 'pad':1,
                   'borderpad':0}

def main(files_dict, ylim = (-60,30), sb_kwds = default_sb_prop):
    # change some legend plot stuff
    plt.rcParams.update(\
        {'legend.frameon':False,
         'legend.numpoints':1,
         'legend.loc':'upper left',
         'legend.fontsize':12})

    # organize info for files / conditions / plot labels into dict

    # pre plot grid
    plt.figure()
    nrow = len(files_dict)
    gs = gridspec.GridSpec(nrow, 1, hspace = 0, wspace = 0,
                           left = 0, right = 1.)

    # plot in loop
    xlims = []
    for i, fdict in enumerate(files_dict):
        plt.subplot(gs[i,0])
        epch = epoch(fdict.pop('abf'), 1, 1)
        epch.set_pading(left = 1000, right = 3000)
        wf = waterfall(epch)
        wf.set_range() # this will change the xlim
        xlims.append(wf._xlim)
        for i, (xs,ys,lvl) in enumerate(wf):
            plt.plot(xs,ys, linewidth = 0.5, **fdict)
            plt.text(min(xs),min(ys), str(lvl), size = 8, va = 'top', ha = 'left')
            if i==0:
                fdict.pop('label')
    
    # post plotting adjust axis so grids are comparable, maybe can do
    # this with some shared axis object? Also hide axis ticks, labels
    # and spines.
    xmin = min([x[0] for x in xlims])
    xmax = max([x[1] for x in xlims])
    for i,ax in enumerate(plt.gcf().axes):
        ax.set_xlim((xmin, xmax)) 
        ax.set_ylim(ylim)
        ax.legend()

        # add a scale bar
        if i==0:
            add_scalebar(ax, bbox_transform = ax.transAxes, **sb_kwds)
        else:
            # hide axis and spines
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)            
            for loc, spine in ax.spines.iteritems():
                    spine.set_color('none')
            
    return plt.gcf()
    
if __name__=='__main__':
    # make list of dicts like so:
    files_dict = [{'abf':abf_reader('2012_08_24_0000.abf'),
                   'color':'black',
                   'label':'control'},
                  {'abf':abf_reader('2012_08_24_0004.abf'),
                   'color':'red',
                   'label':'frf'}]
    fig = main(files_dict)
    fig.set_size_inches(10,8)
    fig.savefig('default_fig_name.png', dpi = 300)
