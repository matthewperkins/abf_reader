import pdb
import os
import numpy as np
from tempfile import mkdtemp

def get_digi_train_period(abr_header):
    digi_trains = np.where(abr_header['ext_epoch_waveform_pulses']\
                 ['nDigitalTrainValue'][0]==4)[0]
    if not np.size(digi_trains)==0:
        return (abr_header['ext_multi-chan_inf']['lEpochPulsePeriod'][0])

def dac_step(level, leng):
    from numpy import tile
    return tile(lvl, leng)

def dac_ramp(rampto, startat, leng):
    from numpy import r_
    from scipy.interpolate import interp1d
    # make interpolation function, prob better way
    return (interp1d([0,leng], [startat,rampto]))

def are_epochs_ramps(abr_header, DAC_num):
    epoch_types =\
        abr_header['ext_epoch_waveform_pulses']['nEpochType'][0][DAC_num]
    return (epoch_types==2)

def find_actv_epchs(abr_header, DAC_num):
    epch_types = abr_header['ext_epoch_waveform_pulses']\
        ['nEpochType'][0, DAC_num]
    actv_epchs = np.nonzero(epch_types)[0]

    # set inactive epchs that are zero length
    eewp = abr_header['ext_epoch_waveform_pulses']
    num_epsds = get_num_episodes(abr_header)
    dur_inits = eewp['lEpochInitDuration'][0]\
        [DAC_num][actv_epchs]
    dur_incrms = eewp['lEpochDurationInc'][0]\
        [DAC_num][actv_epchs]
    NonZDurEpchs = np.where(np.abs(dur_incrms) + np.abs(dur_inits)!=0)
    NonZActvEpchs = actv_epchs[NonZDurEpchs]
    return NonZActvEpchs

def get_num_episodes(abf_header):
    return abf_header['fid_size_info']['lActualEpisodes'][0]

def make_range_array(abf_header, DAC_num):
    num_epsds = get_num_episodes(abf_header)
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    num_epchs = len(actv_epchs)
    range_array = np.tile(np.c_[0:num_epsds],(1,num_epchs))
    return range_array

def get_dp_pad(abf_header):
    SmplsPerEpsd = abf_header['trial_hierarchy']\
        ['lNumSamplesPerEpisode'][0]
    num_chans = get_num_chans(abf_header)
    RowsPerEpsd = (SmplsPerEpsd/num_chans)

    # FINALLY FIGURED OUT HOW THE PRE AND POST HOLDS ARE
    # DETERMINED, THEY ARE EPISODE LENGTH / 64 (MODULO DIV)
    # from the pclamp guide on LTP of all things
    dp_one_side_pad = int(RowsPerEpsd) // int(64)
    return dp_one_side_pad

def make_epch_levels(abf_header, DAC_num):
    eewp = abf_header['ext_epoch_waveform_pulses']
    num_epsds = get_num_episodes(abf_header)
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    num_epchs = len(actv_epchs)

    # construct the levels array
    level_inits = eewp['fEpochInitLevel'][0]\
        [DAC_num][actv_epchs]
    level_incrms = eewp['fEpochLevelInc'][0]\
        [DAC_num][actv_epchs]
    range_array = make_range_array(abf_header, DAC_num)
    tmp_lvl_incrms = np.tile(level_incrms, (num_epsds,1)) * range_array
    levels = np.tile(level_inits,(num_epsds,1))+tmp_lvl_incrms
    return levels

def make_epch_durs(abf_header, DAC_num):
    eewp = abf_header['ext_epoch_waveform_pulses']
    num_epsds = get_num_episodes(abf_header)
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    dur_inits = eewp['lEpochInitDuration'][0]\
        [DAC_num][actv_epchs]
    dur_incrms = eewp['lEpochDurationInc'][0]\
        [DAC_num][actv_epchs]
    range_array = make_range_array(abf_header, DAC_num)
    tmp_dur_incrms = np.tile(dur_incrms, (num_epsds,1)) * range_array
    durs = np.tile(dur_inits,(num_epsds,1))+tmp_dur_incrms
    return (durs)

def get_epch_types(abf_header, DAC_num):
    # epoch types
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    eewp = abf_header['ext_epoch_waveform_pulses']
    epoch_types = eewp['nEpochType'][0]\
        [DAC_num][actv_epchs]
    return (epoch_types)

def get_num_chans(abf_header):
    return np.sum(abf_header['multi-chan_inf']['nADCSamplingSeq'][0]!=-1)
    
def get_total_aquired(abf_header):
    return (abf_header['fid_size_info']['lActualAcqLength'][0])

def sample_rate(header):
    # sample interval in microseconds, so hertz are * 10e6
    sample_rate = \
        1 * 1000000\
        /header['trial_hierarchy']['fADCSampleInterval']\
        /header['trial_hierarchy']['nADCNumChannels']
    sample_rate = sample_rate[0]
    return sample_rate

def get_DAC_len(abf_header):
    tot_aq = get_total_aquired(abf_header)
    num_chns = get_num_chans(abf_header)
    return ( tot_aq // num_chns )

def get_sweep_len(abf_header):
    DAC_ln = get_DAC_len(abf_header)
    num_epsds = get_num_episodes(abf_header)
    return ( DAC_ln // num_epsds )

def get_epsd_len(abf_header):
    swp_ln = sweep_len(abf_header)
    pad = get_dp_pad(abf_header)
    return ( swp_ln - 2*pad )

def make_swp_indxs(abf_header):
    sweep_len = get_sweep_len(abf_header)
    DAC_len = get_DAC_len(abf_header)
    lft = np.r_[0:DAC_len+1:sweep_len][:-1]
    rght = np.r_[0:DAC_len+1:sweep_len][1:]
    sweep_indxs = np.c_[lft,rght]
    return sweep_indxs

def make_epsd_indxs(abf_header):
    sweep_len = get_sweep_len(abf_header)
    sweep_indxs = make_swp_indxs(abf_header)
    epsd_indxs = np.copy(sweep_indxs)
    pad = get_dp_pad(abf_header)
    epsd_indxs[:,0] = sweep_indxs[:,0] + pad
    epsd_len = sweep_len - 2*pad
    epsd_indxs[:,1] = epsd_indxs[:,0]+epsd_len
    return epsd_indxs

def make_epch_indxs(abf_header, DAC_num, nrmd = False):
    #### epoch starts (rel to episode start) ####
    
    # assume left points are fixed in
    # duration changing epochs. I hope this means the maximum duration
    # for an epoch (over the whole trail) must pass before next epoch
    # statrs
    durs = make_epch_durs(abf_header, DAC_num)
    max_dur = np.max(durs,0)
    num_epsds = get_num_episodes(abf_header)
    max_durs = np.tile(max_dur, (num_epsds,1))

    # accumulate durations (add a zero column) to get starts
    strts = np.c_[np.repeat(0,num_epsds), max_durs[:,0:-1]]
    epch_strt_indxs = np.add.accumulate(strts,1)

    #### epoch starts (rel to 0) ####
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    num_epchs = len(actv_epchs)

    epsd_indxs = make_epsd_indxs(abf_header)

    for i_epch in range(num_epchs):
        if nrmd=='sweep':
            epch_strt_indxs[:,i_epch] += (epsd_indxs[0,0])
        elif nrmd=='episode':
            pass
        else:
            epch_strt_indxs[:,i_epch]+=epsd_indxs[:,0]

    epch_end_indxs = epch_strt_indxs+durs

    epch_indx_list = []
    for i_epsd in range(num_epsds):
        tmp_epsd_list=[]
        for i_epch in range(num_epchs):
            t = [epch_strt_indxs[i_epsd,i_epch],
                 epch_end_indxs[i_epsd,i_epch]]
            tmp_epsd_list.append(t)
        epch_indx_list.append(tmp_epsd_list)
    epch_indxs = np.array(epch_indx_list)
    return epch_indxs

class DAC(object):
    def __init__(self, abf_header, DAC_num, cap_pad_ms = 1.5, **kwds):
        from abf_header_dtype import abf_header_dtype
        assert 0 <= DAC_num < 2, "DAC_num must be 0 or 1"
        assert abf_header.dtype is abf_header_dtype, "header has wrong dtype"
        self._num_episodes = get_num_episodes(abf_header)
        self._actv = abf_header['ext_epoch_waveform_pulses']['nWaveformEnable'][0][DAC_num]
        self._epsd_ixs = make_epsd_indxs(abf_header)
        self._epch_ixs = make_epch_indxs(abf_header, DAC_num)
        self._epch_lvls = make_epch_levels(abf_header, DAC_num)
        self._DAC_num = DAC_num
        self._cap_pad_ms = cap_pad_ms
        self._cap_pad_dp = int(cap_pad_ms//1000.*sample_rate(abf_header))
        super(DAC, self).__init__(**kwds)

        self._actv_epchs = find_actv_epchs(abf_header, DAC_num)
        self._num_epchs = len(self._actv_epchs)
        self._cap_pad = False
        self._make_epoch_slices()

    def _make_epoch_slices(self):
        # make epoch slices
        self._epch_slices = []
        for epsd_i in range(self._num_episodes):
            _tmp = []
            for epch_i in range(self._num_epchs):
                _tmp.append(slice(self._epch_ixs[epsd_i, epch_i, 0],
                                  self._epch_ixs[epsd_i, epch_i, 1]))
            self._epch_slices.append(_tmp)

    def pad_cap(self, pad=True):
        if pad is True and self._cap_pad is not True:
            self._epch_ixs[:,:,0]+=self._cap_pad_dp
            self._make_epoch_slices()
            self._cap_pad = True
        elif pad is False and self._cap_pad is True:
            self._epch_ixs[:,:,0]-=self._cap_pad_dp
            self._make_epoch_slices()
            self._cap_pad = False

    def __iter__(self, epoch = 1):
        for epsd_num in range(self._num_episodes):
            yield self._epch_slices[epsd_num][epoch]

class abf_reader(object):
    from mhp_re import yyyy_mm_dd_nnnn
    fnum_re = yyyy_mm_dd_nnnn
    def __init__(self, fname, **kwds):
        from abf_header_dtype import abf_header_dtype
        self._headr_struct = abf_header_dtype
        self.fname = os.path.basename(fname)
        if abf_reader.fnum_re.search(self.fname):
            m = abf_reader.fnum_re.search(self.fname)
            self.fnum = m.groups()[-1]
        else:
            self.fnum = 'NA'
        if os.path.isabs(fname):
            self.path = os.path.dirname(fname)
        else:
            self.path = os.path.dirname(os.path.abspath(fname))
        self.path_file = self.path + os.sep + self.fname
        self.fid = open(self.path_file, 'rb')
        self.read_header()

        # make sure that I have a compatible abf version
        self.verify_version()
        self.addGain()
        self._chan_holder = -1
        self._num_episodes = \
            get_num_episodes(self.header)

        # rewrite the ADC units into a convience variable, trunc to 2 chars
        self._ADC_Units = \
            np.array(self.header['multi-chan_inf']['sADCUnits'][0],
                     dtype = '|S2')
        # rewrite the DAC units into a convience variable, trunc to 2 chars
        self._DAC_Units = \
            np.array(self.header['multi-chan_inf']['sDACChannelUnits'][0],
                     dtype = '|S2')

        # make an atomic size, so that data can be broken up with out
        # segmenting cols(channels)
        if self.header['f_structure']['nDataFormat'][0]==1: #float data
            self.base_size = 4 * self.num_chans() # 4byte size per float
        elif self.header['f_structure']['nDataFormat'][0]==0: #integer data
            self.base_size = 2 * self.num_chans() # 2byte size per int

        # check operation type, continous is 3 episodic stimulation is 5 
        if self.header['fid_size_info']['nOperationMode'][0]==5:
            self.DAC_0 = DAC(self.header, 0)
            self.DAC_1 = DAC(self.header, 1)

        # for covience, make a protocol variable
        self.protocol = self.header['ext_environment_inf']['sProtocolPath'][0].rstrip()

        # to deal with bad
        self.bad_tele = False
        if 'bad_tele' in kwds.keys():
            assert type(kwds['bad_tele']) is bool, 'must be bool'
            self.bad_tele = kwds.pop('bad_tele')
            if self.bad_tele:
                self._fix_bad_tele()

    def verify_version(self):
        FVerNum = self.header['fid_size_info']['fFileVersionNumber']
        ErrMsg = "%s is version %f, this 'prog' only reads abf 1.8" % (self.fname, FVerNum)
        assert (FVerNum>=1.8) & (FVerNum<2.0), ErrMsg

    def hdr_offset(self):
        from abf_header_defs import ABF_BLOCKSIZE
        return ((self.header['f_structure']['lDataSectionPtr'] * ABF_BLOCKSIZE)[0])

    def total_aq(self):
        return (self.header['fid_size_info']['lActualAcqLength'][0])

    def actv_dacs(self):
        return (np.nonzero(self.header['ext_epoch_waveform_pulses']\
            ['nWaveformEnable'][0])[0])
        
    # custom get and set state allow pickle to handel the pickleing of
    # object with out choking on file

    def __getstate__(self):
        odict = self.__dict__.copy() # copy the dict since we change it
        del odict['fid']              # remove filehandle entry
        if 'mm' in odict.keys(): del odict['mm']  # clear memorym when serializing
        return odict

    def __setstate__(self, dict):
        path_file = dict['path_file']
        self.fid = open(path_file, 'rb')

    def read_header(self):
        self.fid.seek(0)
        self.header = np.fromfile(self.fid, self._headr_struct, 1)

    def _read_seq(self):
        return [read for read in self.header['multi-chan_inf']['nADCSamplingSeq'][0] if read != -1]

    def next_chan(self):
        self._chan_holder += 1
        if self._chan_holder > self.num_chans()-1:
            raise StopIteration
        return self._chan_holder

    def num_chans(self):
        return len(self._read_seq())

    def get_chan_name(self, grep_str):
        import re
        from re import compile as recomp
        prog = recomp(grep_str, flags = re.I)
        has_name = []
        self._chan_holder = -1
        try:
            while True:
                chan_indx = self.next_chan()
                chan_name = self.chan_name(chan_indx)
                result = prog.search(chan_name)
                if (result):
                    has_name.append(chan_indx)
        except StopIteration:
            self._chan_holder = -1
            return has_name

    def chan_names(self):
        adc_l = ['adc_' + str(read) for read in np.r_[0:16]]
        chans = self.header['multi-chan_inf']['sADCChannelName'][0]
        sampled_chans = self._read_seq()
        #these list of sampled chans is in the order it was sampled
        for num, sampled_chan in enumerate(sampled_chans):
            print('%-3s' '%-3s' '%-8s' '%-10s' '%-10s' %(num, '-'*3,
                                                         adc_l[sampled_chan],
                                                         '-'*8, chans[sampled_chan]))

    def chan_name(self, chan_no=0):
        chans = self.header['multi-chan_inf']['sADCChannelName'][0]
        sampled_chans = self._read_seq()
        #rstrip removes the trailing white space
        return chans[sampled_chans[chan_no]].rstrip()

    def read_data(self, **kwds):
        '''reads multiplexed data from abfs into an array'''
        ## the times that are asso with discontinuous recording are wonky
        from numpy import float32, int16, memmap

        # prep to establish order of array in memory
        from abf_header_defs import ABF_BLOCKSIZE
        offset = self.hdr_offset()
        total_aq = self.total_aq()
        numchans = self.num_chans()
        ncols = numchans
        nrows = total_aq//numchans

        # handle optional kwds, for subsetting data
        # start_row and needs to be transulated into byte offset
        # other kwds, num_rows or stop_row do not
        if 'start_time' in kwds.keys():
            start_time = kwds.pop('start_time')
            kwds['start_row'] = int(start_time * self.sample_rate())
        if 'stop_time' in kwds.keys():
            stop_time = kwds.pop('stop_time')
            kwds['stop_row'] = int(stop_time * self.sample_rate())
        if 'r_slice' in kwds.keys():
            row_slice = kwds.pop('r_slice')
            start_row = row_slice.start
            offset += (start_row * self.base_size)
            stop_row = row_slice.stop
            # check if start_row is beginning
            if offset!=self.hdr_offset:
                nrows = stop_row - start_row
            elif offset==self.hdr_offset:
                nrows = stop_row
        if 'start_row' in kwds.keys():
            start_row = kwds.pop('start_row')
            offset += (start_row * self.base_size)
        if 'num_rows' in kwds.keys():
            nrows = kwds.pop('num_rows')
        if 'stop_row' in kwds.keys():
            # check if start_row is beginning
            stop_row = kwds.pop('stop_row')
            if offset!=self.hdr_offset:
                nrows = stop_row - start_row
            elif offset==self.hdr_offset:
                nrows = stop_row
         
        #see if is float data
        if self.header['f_structure']['nDataFormat'][0]==1: 
            data = memmap(self.fid, dtype = float32,
                          shape = (nrows,ncols), offset = offset)
            data = np.copy(data)
            return data

        #see if is integer data
        elif self.header['f_structure']['nDataFormat'][0]==0: 
            unscl_data = memmap(self.fid, dtype = int16,
                          shape = (nrows,ncols),
                mode = 'r',offset = offset)
            # now scale data and return
            unscl_data = unscl_data.astype(float32)
            return (self.scale_int_data(unscl_data))

    def scale_int_data(self, data):
        for indx, chan in enumerate(self._read_seq()):
            divis = (self.header['multi-chan_inf']['fInstrumentScaleFactor'][0][chan] * \
                     self.header['multi-chan_inf']['fSignalGain'][0][chan] * \
                     self.header['multi-chan_inf']['fADCProgrammableGain'][0][chan] * \
                     self.addGain[chan])
            mult =  self.header['hardware_inf']['fADCRange'][0] \
                   / self.header['hardware_inf']['lADCResolution'][0]
            offs = self.header['multi-chan_inf']['fInstrumentOffset'][0][chan] - \
                   self.header['multi-chan_inf']['fSignalOffset'][0][chan]
            data[:,indx] = data[:,indx] / divis * mult + offs
        return data

    def addGain(self):
        '''method helps with scaling'''
        self.addGain = self.header['ext_environment_inf']['nTelegraphEnable'][0] * \
            self.header['ext_environment_inf']['fTelegraphAdditGain'][0]
        self.addGain = np.where(self.addGain==0, 1, self.addGain)

    def _fix_bad_tele(self):
        '''hack to work around abf 1.8s created with clampfit10'''
        self.fid_w = open(self.path_file, 'r+b')
        self.fid_w.seek(0)
        # just disable the telegraph gains by setting to zero
        self.header['ext_environment_inf']['nTelegraphEnable'][0]*=0
        self.fid_w.write(self.header)
        self.fid_w.close()
        del self.fid_w

    def get_synch_array(self):
        from abf_header_defs import ABF_BLOCKSIZE
        self.fid.seek(self.header['f_structure']['lSynchArrayPtr'][0] * ABF_BLOCKSIZE)
        synch_array_dtype = [('start', np.int32), ('length', np.int32)]
        synch_array = np.fromfile(self.fid,
                                  synch_array_dtype,
                                  self.header['f_structure']['lSynchArraySize'])
        sl = []
        for strt, length in synch_array:
            sl.append([strt, length])
        sa = np.array(sl)
        return sa

    def sample_rate(self):
        try:
            self._sample_rate
            return self._sample_rate
        except AttributeError:
          # sample interval in microseconds, so hertz are * 10e6
          self._sample_rate = \
          1 * 1000000\
          /self.header['trial_hierarchy']['fADCSampleInterval']\
          /self.header['trial_hierarchy']['nADCNumChannels']
          self._sample_rate = self._sample_rate[0]
          return self._sample_rate

    def xstep(self):
        try:
            self._xstep
            return self._xstep
        except AttributeError:
          # sample interval in microseconds, so hertz are * 10e6
          return 1/self.sample_rate()

    def start_time(self, new_time = None):
        if new_time is None:
            try:
                return self._file_start_time
            except AttributeError:
                from datetime import datetime
                self._File_Time = {}
                yyyymmdd = str(self.header['fid_size_info']['lFileStartDate'][0])
                self._File_Time['year'] = int(yyyymmdd[0:4])
                self._File_Time['month'] = int(yyyymmdd[4:6])
                self._File_Time['day'] = int(yyyymmdd[-2:])

                # 'lFileStartTime is in seconds. do some division for time
                seconds_time = self.header['fid_size_info']['lFileStartTime'][0]
                self._File_Time['hour'] = seconds_time//(60*60)
                self._File_Time['minute'] = (seconds_time%(60*60))//60
                self._File_Time['second'] = (seconds_time%(60*60))%60
                self._File_Time['microsecond'] = \
                  int(self.header['environment_inf']['nFileStartMillisecs'][0]\
                  * 1000)

                #for reading self._File_Time = t_d
                t_d = self._File_Time
                self._file_start_time = datetime(t_d['year'],\
                                          t_d['month'] , t_d['day'],\
                                          t_d['hour'],t_d['minute'],\
                                          t_d['second'],t_d['microsecond'])
                return self._file_start_time
        else:
            from datetime import datetime
            from math import floor
            assert type(new_time) is datetime
            NewTimeDayStart = datetime(new_time.year,
                                 new_time.month,
                                 new_time.day)
            FileStartDate = np.zeros(1, dtype = np.int32)
            FileStartDate[:] = new_time.year * 10**4 + \
                new_time.month * 10**2 +\
                new_time.day
            FileStartTime = np.zeros(1,np.int32)
            Seconds = (new_time - NewTimeDayStart).total_seconds()
            FileStartTime[:] = floor(Seconds)
            FileStartMilli = np.zeros(1,np.int32)
            FileStartMilli[:] = int((Seconds - floor(Seconds))*1000)
            self.header['fid_size_info']['lFileStartDate'][0] = FileStartDate
            self.header['fid_size_info']['lFileStartTime'][0] = FileStartTime
            self.header['environment_inf']['nFileStartMillisecs'][0] = FileStartMilli
            self.fid_w = open(self.path_file, 'r+b')
            self.fid_w.seek(0)
            self.fid_w.write(self.header)
            self.fid_w.close()
            del self.fid_w
            try:
                del self._file_start_time
            except AttributeError:
                pass

    def stop_watch_time(self):
        return int(self.header['fid_size_info']['lStopwatchTime'][0])

if __name__=='__main__':
    import os
    pth = os.path.join(os.environ.get("LABDIR"),
                       'B44_B48_experiments','cng_current',
                       'cbi2_current_time_course','expts',
                       '2013_12_10','clmpx_bin')
    abr = abf_reader(os.path.join(pth,'2013_12_10_0023.abf'))
    d = abr.DAC_0


