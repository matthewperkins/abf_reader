from abf_reader import *
import matplotlib.pyplot as plt
import warnings
import numpy as np

class epoch(object):
    def __init__(self, abf, DAC_num, epoch_num, **kwds):
        self.abf = abf
        # convience variable for sample rate
        self.sr = sample_rate(self.abf.header)
        self.DAC_num = DAC_num
        self.epoch_num = epoch_num
        self._levels = make_epch_levels(abf.header, DAC_num)[:,epoch_num]
        self.epch_indxs = make_epch_indxs(abf.header,
                                     DAC_num,
                                     nrmd = 'sweep')[:,epoch_num]
        self.set_pading(left = 0, right = 0)
        self._num_epchs = len(self.epch_indxs)
        # need to have sweep indexs, bc epoch indexs are normalized to these
        self._swp_indxs = make_swp_indxs(self.abf.header)
        self._num_swps = get_num_episodes(self.abf.header)
        self.__left_pad = 0
        self.__right_pad = 0
        super(epoch, self).__init__(**kwds)

    def set_pading(self, left = 0, right = 0):
        ''' can change how much of the adjacent epochs to include when
        returning data, Should be able to take negative arguments for
        left and right, if you want to cut into epoch. Bounds
        checking here on negative values not set up '''
        # set up arrays for padding 
        # padding is limited to the extent of the sweep
        sweep_len = get_sweep_len(self.abf.header)
        self.padded_epch_indxs = np.copy(self.epch_indxs)
        # left indexs, pad than zero out negatives
        self.padded_epch_indxs[:,0] -= left
        self.padded_epch_indxs[:,0][np.where(self.padded_epch_indxs[:,0]<0)[0]]=0
        self.padded_epch_indxs[:,0][np.where(self.padded_epch_indxs[:,0]>sweep_len)]=sweep_len
        # right indxes, pad than set overlong to sweep len
        self.padded_epch_indxs[:,1] += right
        self.padded_epch_indxs[:,1][np.where(self.padded_epch_indxs[:,1]>sweep_len)]=sweep_len
        self.padded_epch_indxs[:,1][np.where(self.padded_epch_indxs[:,1]<0)[0]]=0
        # now loop through just to check
        for i, (strt, stop) in enumerate(self.padded_epch_indxs):
            assert strt!=stop, "Epoch Padding is BAD: Start Idx (%d) == Stop Idx (%d) at Sweep %d" % (strt, stop, i)

    def _data(self, swp_idx):
        strt = self._swp_indxs[swp_idx,0] +\
                    self.padded_epch_indxs[swp_idx,0]
        stop = self._swp_indxs[swp_idx,0] +\
                    self.padded_epch_indxs[swp_idx,1]
        if stop<strt:
            # if stop is before start, read and reverse
            warn_msg = "\n\nStart Epoch Idx (%d) is after Stop Epoch Idx (%d), for sweep (%d), data are time inverted\n\n" % (strt, stop, swp_idx)
            warnings.warn(warn_msg)
            return (self.abf.read_data(start_row = stop, stop_row = strt)[::-1])
        elif strt<stop:
            return (self.abf.read_data(start_row = strt, stop_row = stop))

    def _lvl(self, swp_idx):
        return self._levels[swp_idx]
        
    def __iter__(self):
        return self.gen_iter()

    def gen_iter(self):
        # using generator
        for swp_idx in range(self._num_swps):
            yield self._data(swp_idx)

class waterfall(object):
    def __init__(self, abf_epoch,
                 channo = 0, xprcnt = 110, xoffset_sec = False,
                 yoffset = 0, selectd = False,
                 **kwds):
        self.abf_epoch = abf_epoch
        self.channo = channo
        self.xprcnt = xprcnt
        self.yoffset = yoffset
        self._num_swps = self.abf_epoch._num_swps
        self.swp_len_dp = len(abf_epoch.__iter__().next()[:,self.channo])
        self.sr = sample_rate(abf_epoch.abf.header)
        self.swp_len_sec = self.swp_len_dp/self.sr
        if xprcnt is not False:
            self.offset = (self.swp_len_sec) * xprcnt/100.
            self.interstitial = self.offset - self.swp_len_sec
            self._fixd_offset = False
        if xoffset_sec is not False:
            # sec an x offset in seconds
            self.offset = xoffset_sec
            self.interstitial = 0
            self._fixd_offset = True
        self.xs = np.linspace(0, self.swp_len_sec, self.swp_len_dp)
        if selectd is False:
            self._swps = np.arange(self._num_swps)
        else:
            assert issubclass(type(selectd), np.ndarray), "selectd is type %s, should be np.ndarray" % type(selectd).__name__
            assert len(selectd.shape)==1
            self._swps = selectd
        # just to not break things, too bad that is 
        self.set_range = self.x_range
        super(waterfall, self).__init__(**kwds)

    def __iter__(self, **kwds):
        return self.gen_iter(**kwds)

    def x_range(self):
        if self._fixd_offset:
            self._xlim = (\
                0, len(self._swps) * self.offset + self.swp_len_sec - self.offset)
        else:
            self._xlim =\
                (0, len(self._swps) * self.swp_len_sec +\
                     (self.interstitial * self._num_iter_sweeps))
        return self._xlim

    def gen_iter(self):
        for i, swp_idx in enumerate(self._swps):
            xs = self.xs + self.offset*i
            ys = self.abf_epoch._data(swp_idx)[:,self.channo] + self.yoffset*i
            lvl = self.abf_epoch._lvl(swp_idx)
            yield xs, ys, lvl
        self.xlim = (0, max(xs))

if __name__=='__main__':
    ''' debug script '''
    import os
    labdir = os.environ.get("LABDIR")
    abfpath = os.path.join(labdir,
                 'B44_B48_experiments',
                 'b48_fI_summary',
                 'frf', '2012_08_24_0004.abf')
    abf = abf_reader(abfpath)
    epch = epoch(abf, 1, 1)
    epch.set_pading(left = 1000, right = 1000)
    for xs,ys,lvl in waterfall(epch, selectd = np.r_[0:18:2]):
        plt.plot(xs,ys, '-k', linewidth = 0.5)
    plt.show()
