import os
from abf_reader import *
import matplotlib.pyplot as plt
from matt_axes_cust import xscale, yscale, clean_axes, no_spines
import pdb

########## load fixture data ##########
data_path = os.path.join(os.environ.get('LABDIR'),
                         'B44_B48_experiments','b48_fI_summary','control')
fI = abf_reader(os.path.join(data_path, '2012_09_08_0009.abf'))
########## end load fixture data ##########

#################### time steps and DACs #################### 
DAC_num = 1                       
ADC_num = 0
xstep = 1/sample_rate(fI.header)
sr = sample_rate(fI.header)
#################### end time steps and DACs ################

#################### PLOT fI ####################
for swp_num, swp in enumerate(make_swp_indxs(fI.header)):
    swp_len = swp[1]-swp[0]
    swp_x_offset = 1.4 * swp_num
    swp_y_offset = 0 * swp_num
    xs = np.r_[0:swp_len/sr:xstep] + swp_x_offset
    ys = fI.read_data()[slice(swp[0],swp[1]),0] + swp_y_offset
    # make sure xs arnt too long, with float fups
    xs = xs[0:len(ys)]
    plt.plot(xs,ys, linewidth = 0.3, color = 'black')
    xmax = np.max(xs)
ax = plt.gca()
fig =plt.gcf()
ax.set_ylim((-65,25))
ax.set_xlim((0,xmax))

clean_axes(ax)
no_spines(ax)
# from matt_scalebars import mk_scale_bar
# mk_scale_bar(ax, horz_length_data = 1, vert_length_data = 20,
#              vert_units = 'mV', horz_units = 'sec', arrangement = 'qoph',
#              axescoords = (0.05,0.8), linewidth = 3, vpad = -0.03)

from matt_axes_cust import xscale, yscale
xscale(ax, 2, xloc_axs = 0.05, xalgn = 'left')
yscale(ax, 20, xloc_axs = 0.05, xalgn = 'left')
fig.set_size_inches((11,3))
fig.savefig('nice_fI.png', dpi = 110)
#################### END PLOT fI ####################

#################### PLOT latn I ####################
actv_epch = 1
pre_ms = 50
post_ms = 200
all_epch_indxs = make_epch_indxs(fI.header, DAC_num)
for swp_num, epch_indxs in enumerate(all_epch_indxs):
    plt_ln = pre_ms + post_ms
    swp_x_offset = 0 * swp_num
    swp_y_offset = 3 * swp_num
    xs = np.r_[0:plt_ln:xstep*1000] + swp_x_offset
    start = epch_indxs[actv_epch][0] - pre_ms/xstep/1000
    end = epch_indxs[actv_epch][0] + post_ms/xstep/1000
    ys = fI.read_data()[slice(start,end),ADC_num] + swp_y_offset
    # make sure xs arnt too long, with float fups
    xs = xs[0:len(ys)]
    plt.plot(xs,ys, linewidth = 0.5, color = 'black', zorder = swp_num-len(all_epch_indxs))
    xmax = np.max(xs)
ax = plt.gca()
fig =plt.gcf()
ax.set_ylim((-65,25+swp_y_offset*0.8))
ax.set_xlim((0,xmax))

clean_axes(ax)
no_spines(ax)
# from matt_scalebars import mk_scale_bar
# mk_scale_bar(ax, horz_length_data = 1, vert_length_data = 20,
#              vert_units = 'mV', horz_units = 'sec', arrangement = 'qoph',
#              axescoords = (0.05,0.8), linewidth = 3, vpad = -0.03)

from matt_axes_cust import xscale, yscale
xscale(ax, 50, xloc_axs = 0.05, xalgn = 'left', unitstr = 'msec')
yscale(ax, 20, xloc_axs = 0.05, xalgn = 'left')
fig.set_size_inches((11,5))
fig.savefig('nice_fLatn.png', dpi = 150)

