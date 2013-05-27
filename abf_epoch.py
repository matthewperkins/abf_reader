from abf_reader import *
import matplotlib.pyplot as plt

class epoch(object):
    def __init__(self, abf, DAC_num, epoch_num, **kwds):
        self.abf = abf
        self.DAC_num = DAC_num
        self.epoch_num = epoch_num
        self.levels = make_epch_levels(abf.header, DAC_num)
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
        # right indxes, pad than set overlong to sweep len
        self.padded_epch_indxs[:,1] += right
        self.padded_epch_indxs[:,1][np.where(self.padded_epch_indxs[:,1]>sweep_len)]=sweep_len

    def data(self, swp_idx):
        strt = self._swp_indxs[swp_idx,0] +\
                    self.padded_epch_indxs[swp_idx,0]
        stop = self._swp_indxs[swp_idx,0] +\
                    self.padded_epch_indxs[swp_idx,1]
        return (self.abf.read_data(start_row = strt, stop_row = stop))
        
    def __iter__(self):
        return self.gen_iter()

    def gen_iter(self):
        # using generator
        for swp_idx in range(self._num_swps):
            yield self.data(swp_idx)

class waterfall(object):
    def __init__(self, abf_epoch,
                 channo = 0, xprcnt = 110, yoffset = 0, **kwds):
        self.abf_epoch = abf_epoch
        self.channo = channo
        self.xprcnt = xprcnt
        self.yoffset = yoffset
        self._num_swps = self.abf_epoch._num_swps
        self.swp_len_dp = len(abf_epoch.__iter__().next()[:,self.channo])
        self.sr = sample_rate(abf_epoch.abf.header)
        self.swp_len_sec = self.swp_len_dp/self.sr
        self.offset = (self.swp_len_sec) * xprcnt/100.
        self.interstitial = self.offset - self.swp_len_sec
        self.xs = np.linspace(0, self.swp_len_sec, self.swp_len_dp)
        super(waterfall, self).__init__(**kwds)

    def __iter__(self, **kwds):
        return self.gen_iter(**kwds)

    def set_range(self, start = 0, stop = -1, step = 1):
        stop = stop if stop>0 else self._num_swps
        self._stop = stop
        self._start = start
        self._step = step
        self._num_iter_sweeps =\
            len(range(self._start, self._stop, self._step))
        self._xlim =\
            (0, self._num_iter_sweeps * self.swp_len_sec +\
                 (self.interstitial * self._num_iter_sweeps))

    def gen_iter(self):
        for i, swp_idx in enumerate(range(self._start,
                                          self._stop,
                                          self._step)):
            xs = self.xs + self.offset*i
            ys = self.abf_epoch.data(swp_idx)[:,self.channo]
            yield xs, ys
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
    for xs,ys in waterfall(epch):
        plt.plot(xs,ys, '-k', linewidth = 0.5)
    plt.show()
