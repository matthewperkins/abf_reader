import matplotlib.gridspec as gridspec
from matt_axes_cust import no_spines
import pdb
import numpy as np
import matplotlib.pyplot as plt
from signal_tools import bwfiltfilt
from mp_scale_bars import AnchoredScaleBar, DraggableScaleBar
import os
import copy

def trace_axes(ax):
    from matplotlib.pyplot import setp
    for loc, spine in ax.spines.items():
        spine.set_color('none') # don't draw spine
    setp(ax.get_xticklabels(), visible = False)
    setp(ax.get_yticklabels(), visible = False)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

class abf_chunker_plotter(object):
    def __init__(self, _abf_chunker):
        self.num_chunks = _abf_chunker.num_chunk_rows
        self._abf_chunker = _abf_chunker
        self.set_fig_size_inches((8,3))
        self._ylim = {}
        self.set_linewidth()

    def set_ylim(self, **kwds):
        for kwd in kwds.keys():
            if kwd in ['cells','neurgrms']:
                self._ylim[kwd] = kwds[kwd]

    def set_fig_size_inches(self, size):
        #size in inches / number chunks(float)
        self.indiv_fig_width = float(size[0]) /\
                                 (self._abf_chunker.frct_chunks)
        self.indiv_fig_height = float(size[1])

    def set_linewidth(self, width = False):
        ''' trail and error, seems that narrower pens can be used for shorter
        segments, and thick pen must be used for longer segments'''
        if width!=False:
            self.linethickness = width
            return
        
        total_time = self._abf_chunker.get_range_seconds()
        sample_rate = self._abf_chunker.abr.sample_rate()
        wickets = np.r_[0:110:10, 120:240:20, 220:1050:50]
        thickness = np.r_[0.4:0.8:0.005]
        left_wicket = np.where(wickets<total_time)[0][-1]
        self.linethickness = thickness[left_wicket-1]
        
    def plot(self, dpi, cell_list, neurgrm_list, **kwds):

        # sample scale pars in keyword
        # kwds['scale_bars'] = {'shared':{'bbox_to_anchor' : (0.3,0.4),
        #                                 'xsep' : 1,
        #                                 'ysep' : 1,
        #                                 'pad' : 1,
        #                                 'prop' : FP(size=6),
        #                                 'borderpad':0},
        #                       'cell':[{'sizex' : 0.5,
        #                                'labelx' : '0.5 sec',
        #                                'sizey' : 20,
        #                                'labely' : '20 nA'}]*num_cells,
        #                       'neurgrm':[{'sizex' : 0.5,
        #                                   'labelx' : '0.5 sec',
        #                                   'sizey' : 20,
        #                                   'labely' : '20 nA'}]*num_neurgrms}
        

        # most of the plot will be for intracell

        # the ratio that cell signals will be bigger than neurogram signals
        if 'cell_bias' in kwds.keys():
            cell_bias = kwds['cell_bias']
        else:
            cell_bias = 4/1.

        # calc the fraction of the fig that will be for cells vs neurgrms,
        # using above fraction
        nrgrm_height = len(neurgrm_list)
        cell_height = len(cell_list) * cell_bias 
        tot_height =  nrgrm_height + cell_height
        
        celltop = 1
        cellbtm = 1 - (cell_height / tot_height)

        #neurograms will only take up about bottom thrid
        nrgtop = cellbtm - 0.025 # give a small (vertical 2.5% pad)
        nrgbtm = 0

        #if there is a image_type keyword, pop and set
        if 'image_type' in kwds.keys():
            image_type = kwds.pop('image_type')
        else:
            image_type = 'png'

        #if there is a image_name keyword, pop and set
        if 'image_name' in kwds.keys():
            image_name = kwds.pop('image_name')
        else:
            image_name = 'tmp'

        # if there is a color keyword, pop and set
        if 'color' in kwds.keys():
            col = kwds.pop('color')
        else:
            col = 'black'

        # if there is a color keyword, pop and set
        if 'plt_kwargs' not in kwds.keys():
            kwds['plt_kwargs'] = {'color':col,
                                  'linewidth':self.linethickness}
        
        import subprocess
        tmp_files = []
        for i, d in enumerate(self._abf_chunker):
            fig = plt.figure()
            width_frct = float(len(d))/self._abf_chunker.chunk_row_size
            gs_cells = gridspec.GridSpec(len(cell_list), 1)
            gs_cells.update(left = 0, right = 1,
                            top = celltop, bottom = cellbtm, hspace = 0.05)
            gs_neurgrm = gridspec.GridSpec(len(neurgrm_list), 1)
            gs_neurgrm.update(left = 0, right = 1,
                              top = nrgtop, bottom = nrgbtm, hspace = 0.05)
            if 'scale_bars' in kwds.keys():
                # make a copy to avoid shitty side effects
                sb = copy.deepcopy(kwds['scale_bars'])
                dsbs = []
                # need to adjust the x to data points here
                try: 
                    for SclBarDict in sb['cell']:
                        try:
                            assert SclBarDict['correctdx'] is True
                        except (KeyError, AssertionError):
                            SclBarDict['sizex'] *= self._abf_chunker.abr.sample_rate()
                            SclBarDict['xstep'] = 1/self._abf_chunker.abr.sample_rate()
                            SclBarDict['correctdx'] = True
                    # kludge need to pop this temp keyword off again                            
                    for SclBarDict in sb['cell']:
                        try:
                            SclBarDict.pop('correctdx')
                        except KeyError:
                            pass
                except KeyError:
                    pass
                try:
                    for SclBarDict in sb['neurgrm']:
                        try:
                            assert SclBarDict['correctdx'] is True
                        except (KeyError, AssertionError):
                            SclBarDict['sizex'] *= self._abf_chunker.abr.sample_rate()
                            SclBarDict['xstep'] = 1/self._abf_chunker.abr.sample_rate()
                            SclBarDict['correctdx'] = True
                    # kludge need to pop this temp keyword off again                            
                    for SclBarDict in sb['neurgrm']:
                        try:
                            SclBarDict.pop('correctdx')
                        except KeyError:
                            pass
                except KeyError:
                    pass
            
            for cell_num, cell in enumerate(cell_list):
                data = d[:,cell]
                ax = plt.subplot(gs_cells[cell_num, 0])
                if 'axisbg' in kwds.keys():
                    ax.set_axis_bgcolor(kwds['axisbg'])
                    trace_axes(ax)
                else:
                    ax.set_axis_off()
                if 'lp_filt' in kwds.keys():
                    if kwds['lp_filt'][cell_num] is False:
                        pass
                    else:
                        data = bwfiltfilt(data,
                                          self._abf_chunker.abr.sample_rate(),
                                          kwds['lp_filt'][cell_num])
                    if np.any(np.isnan(data)):
                        raise ValueError('data has NaNs, probably filter is bad')
                if 'scale_bars' not in kwds.keys():
                    if type(kwds['plt_kwargs']) is list:
                        ax.plot(data, **kwds['plt_kwargs'][cell_num])
                    else:
                        ax.plot(data, **kwds['plt_kwargs'])
                ax.set_ylim(self._ylim['cells'][cell_num])
                ax.set_xlim((0,len(d)))
                if 'scale_bars' in kwds.keys():
                    ASB = AnchoredScaleBar(plt.gca().transData,
                                           plt.gcf().dpi_scale_trans,
                                           bbox_transform = plt.gca().transAxes,
                                           **dict(sb['shared'], **sb['cell'][cell_num]))
                    ax.add_artist(ASB)
            for neurgrm_num, neurgrm in enumerate(neurgrm_list):
                data = d[:,neurgrm]
                ax = plt.subplot(gs_neurgrm[neurgrm_num, 0])
                if 'axisbg' in kwds.keys():
                    ax.set_axis_bgcolor(kwds['axisbg'])
                    trace_axes(ax)
                else:
                    ax.set_axis_off()
                if 'scale_bars' not in kwds.keys():
                    ax.plot(data, **kwds['plt_kwargs'])
                ax.set_ylim(self._ylim['neurgrms'][neurgrm_num])
                ax.set_xlim((0,len(d)))
                if 'scale_bars' in kwds.keys():
                    ASB = AnchoredScaleBar(plt.gca().transData,
                                           plt.gcf().dpi_scale_trans,
                                           bbox_transform = plt.gca().transAxes,
                                           **dict(sb['shared'].items(),**sb['neurgrm'][neurgrm_num]))
                    ax.add_artist(ASB)
            print(fig.get_size_inches())
            del d
            filname = "%s_%03d.%s" % ('tmp', i, image_type)
            print(filname)
            if 'scale_bars' in kwds.keys():
                scale_bar = kwds.pop('scale_bars')
                plt.gcf().set_size_inches(( self.indiv_fig_width * width_frct,
                                            self.indiv_fig_height))
                plt.gcf().savefig(filname, dpi = dpi, transparent = True)
                plt.close()
            else:
                fig.set_size_inches(( self.indiv_fig_width * width_frct, self.indiv_fig_height))
                if 'axisbg' in kwds.keys():
                    fig.patch.set_color(kwds['axisbg'])
                    fig.savefig(filname, dpi = dpi,
                                facecolor = kwds['axisbg'],
                                edgecolor = kwds['axisbg'])
                else:
                    fig.savefig(filname, dpi = dpi, transparent = True)
            ### the order that these are cleared is very important, otherwise have memleaks.
            tmp_files.append(filname)
            plt.clf()
            del ax
            plt.close(fig)
            del fig
            
        if self._abf_chunker.num_chunk_rows<1:
            image_name = "%s.%s" % (image_name, image_type)
            if os.path.exists(image_name):
                os.remove(image_name)
            os.rename(filname, image_name)
            return 0

        # this is so cool, calling ImageMagic as a subprocess!!

        # the montage command should be same more or less always
        # montage -tile x1 -background transparent -mode Concatenate `ls *.png` out.png
        montage_cmmnd = ['gm', 'montage','-tile','x1','-background','transparent','-mode','Concatenate']
            
        # add the files to stitch
        montage_cmmnd.extend(tmp_files)
        # add the output file name
        image_name = "%s.%s" % (image_name, image_type)
        print(image_name)
        montage_cmmnd.extend([image_name])

        # do the command with a subprocess
        subprocess.call(montage_cmmnd)


        # remove the temporary files
        tmp_files = map(os.path.abspath, tmp_files)
        map(os.remove, tmp_files)

        # display
        # make the display command platform specific
        # import platform
        # out_path = os.path.abspath('out.png')
        # if platform.system() == 'Windows':
        #     dsply_cmnd = ['gmdisplay', out_path]
        # else:
        #     dsply_cmnd = ['display', out_path]
        # subprocess.call(dsply_cmnd)
        # return os.path.abspath('out.png')
