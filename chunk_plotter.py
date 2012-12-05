import matplotlib.gridspec as gridspec
from matt_axes_cust import no_spines
import pdb
import numpy as np
import matplotlib.pyplot as plt
import os

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

        # most of the plot will be for intracell

        # the ratio that cell signals will be bigger than neurogram signals
        cell_bias = 2/1.

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

        #if there is a filetype keyword, pop and set
        if 'filetype' in kwds.keys():
            filetype = kwds.pop('filetype')
        else:
            filetype = 'png'
        
        import subprocess
        tmp_files = []
        for i, d in enumerate(self._abf_chunker):
            fig = plt.figure()
            fig.patch.set_alpha(0)
            width_frct = float(len(d))/self._abf_chunker.chunk_row_size
            fig.set_size_inches(( self.indiv_fig_width * width_frct, self.indiv_fig_height))
            gs_cells = gridspec.GridSpec(len(cell_list), 1)
            gs_cells.update(left = 0, right = 1,
                            top = celltop, bottom = cellbtm, hspace = 0.05)
            gs_neurgrm = gridspec.GridSpec(len(neurgrm_list), 1)
            gs_neurgrm.update(left = 0, right = 1,
                              top = nrgtop, bottom = nrgbtm, hspace = 0.05)
            for cell_num, cell in enumerate(cell_list):
                data = d[:,cell]
                ax = plt.subplot(gs_cells[cell_num, 0])
                if 'yscale' in kwds.keys():
                    ybar = kwds.pop('yscale')
                    print('ybar')
                    print(ybar)
                    ax.plot([0,0],[0,ybar], color = 'red', linewidth = 1)
                if 'xscale' in kwds.keys():
                    xbar = kwds.pop('xscale')*\
                        self._abf_chunker.abr.sample_rate()
                    print('xbar')
                    print(xbar)
                    ax.plot([0,xbar],[0,0], color = 'red', linewidth = 1)

                ax.set_axis_off()
                ax.plot(data, linewidth = self.linethickness, color = 'black')
                ax.set_ylim(self._ylim['cells'][cell_num])
                ax.set_xlim((0,len(d)))
            for neurgrm_num, neurgrm in enumerate(neurgrm_list):
                data = d[:,neurgrm]
                ax = plt.subplot(gs_neurgrm[neurgrm_num, 0])
                ax.set_axis_off()
                ax.plot(data, linewidth = self.linethickness, color = 'black')
                ax.set_ylim(self._ylim['neurgrms'][neurgrm_num])
                ax.set_xlim((0,len(d)))
            print(fig.get_size_inches())
            del d
            filname = "%s_%03d.%s" % ('tmp', i, filetype)
            print(filname)
            tmp_files.append(filname)
            fig.savefig(filname, dpi = dpi)
            ### the order that these are cleared is very important, otherwise have memleaks.
            plt.clf()
            del ax
            plt.close(fig)
            del fig

        # this is so cool, calling ImageMagic as a subprocess!!

        # the montage command should be same more or less always
        # montage -tile x1 -background transparent -mode Concatenate `ls *.png` out.png
        montage_cmmnd = ['montage','-tile','x1','-background','transparent','-mode','Concatenate']
            
        # add the files to stitch
        montage_cmmnd.extend(tmp_files)
        # add the output file name
        if 'image_name' in kwds.keys():
            image_name = kwds.pop('image_name')
        else:
            image_name = 'out.png'
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
