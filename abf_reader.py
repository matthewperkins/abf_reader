import pdb
import os
import numpy as np
from tempfile import mkdtemp

def dac_step(level, leng):
    from numpy import tile
    return tile(lvl, leng)

def dac_ramp(rampto, startat, leng):
    from numpy import r_
    from scipy.interpolate import interp1d
    # make interpolation function, prob better way
    return (interp1d([0,leng], [startat,rampto]))

def find_actv_epchs(abr_header, DAC_num):
    epch_types = abr_header['ext_epoch_waveform_pulses']\
        ['nEpochType'][0, DAC_num]
    return np.nonzero(epch_types)[0]

def get_num_episodes(abf_header):
    return abf_header['trial_hierarchy']['lEpisodesPerRun'][0]

def make_range_array(abf_header, DAC_num):
    num_epsds = abf_header['trial_hierarchy']['lEpisodesPerRun'][0]
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    num_epchs = len(actv_epchs)
    range_array = np.transpose(np.tile(np.c_[0:num_epsds],
                                       (num_epchs,1)))
    return range_array

def get_dp_pad(abf_header):
    SmplsPerEpsd = abf_header['trial_hierarchy']\
        ['lNumSamplesPerEpisode'][0]
    num_chans = get_num_chans(abf_header)
    RowsPerEpsd = (SmplsPerEpsd/num_chans)

    # FINALLY FIGURED OUT HOW THE PRE AND POST HOLDS ARE
    # DETERMINED, THEY ARE EPISODE LENGTH / 64 (MODULO DIV)
    # from the pclamp guide on LTP of all things
    dp_one_side_pad = int(RowsPerEpsd) / int(64)
    return dp_one_side_pad

def make_epch_levels(abf_header, DAC_num):
    eewp = abf_header['ext_epoch_waveform_pulses']
    num_epsds = abf_header['trial_hierarchy']['lEpisodesPerRun'][0]
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
    num_epsds = abf_header['trial_hierarchy']['lEpisodesPerRun'][0]
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
    read_chans =  filter(lambda read: \
        read != -1, abf_header['multi-chan_inf']['nADCSamplingSeq'][0])
    return len(read_chans)

def get_total_aquired(abf_header):
    return (abf_header['fid_size_info']['lActualAcqLength'][0])

def get_DAC_len(abf_header):
    tot_aq = get_total_aquired(abf_header)
    num_chns = get_num_chans(abf_header)
    return ( tot_aq / num_chns )

def get_sweep_len(abf_header):
    DAC_ln = get_DAC_len(abf_header)
    num_epsds = get_num_episodes(abf_header)
    return ( DAC_ln / num_epsds )

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

def make_epch_indxs(abf_header, DAC_num):
    #### epoch starts (rel to episode start) ####
    
    # assume left points are fixed in
    # duration changing epochs. I hope this means the maximum duration
    # for an epoch (over the whole trail) must pass before next epoch
    # statrs
    durs = make_epch_durs(abf_header, DAC_num)
    max_dur = np.max(durs,0)
    num_epsds = abf_header['trial_hierarchy']['lEpisodesPerRun'][0]
    max_durs = np.tile(max_dur, (num_epsds,1))

    # accumulate durations (add a zero column) to get starts
    strts = np.c_[np.repeat(0,num_epsds), max_durs[:,0:-1]]
    epch_strt_indxs = np.add.accumulate(strts,1)

    #### epoch starts (rel to 0) ####
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    num_epchs = len(actv_epchs)

    epsd_indxs = make_epsd_indxs(abf_header)
    
    for i_epch in range(num_epchs):
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

class abf_reader(object):
    def __init__(self, fname):
        from abf_header_dtype import abf_header_dtype
        self._headr_struct = abf_header_dtype
        self.fname = os.path.basename(fname)
        if os.path.isabs(fname):
            self.path = os.path.dirname(fname)
        else:
            self.path = os.path.dirname(os.path.abspath(fname))
        self.path_file = self.path + os.sep + self.fname
        self.fid = file(self.path_file, 'rb')
        self.read_header()
        self.addGain()
        self._chan_holder = -1
        self._num_episodes = \
            self.header['trial_hierarchy']['lEpisodesPerRun'][0]

        # rewrite the ADC units into a convience variable, trunc to 2 chars
        self._ADC_Units = \
            np.array(self.header['multi-chan_inf']['sADCUnits'][0],
                     dtype = '|S2')
        # make an atomic size, so that data can be broken up with out
        # segmenting cols(channels)
        if self.header['f_structure']['nDataFormat'][0]==1: #float data
            self.base_size = 4 * self.num_chans() # 4byte size per float
        elif self.header['f_structure']['nDataFormat'][0]==0: #integer data
            self.base_size = 2 * self.num_chans() # 2byte size per int

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
        self.fid = file(path_file, 'rb')

    def read_header(self):
        self.fid.seek(0)
        self.header = np.fromfile(self.fid, self._headr_struct, 1)

    def _read_seq(self):
        return filter(lambda read: \
            read != -1, self.header['multi-chan_inf']['nADCSamplingSeq'][0])

    def next_chan(self):
        self._chan_holder += 1
        if self._chan_holder > self.num_chans()-1:
            raise StopIteration
        return self._chan_holder

    def num_chans(self):
        return len(self._read_seq())

    def get_chan_name(self, grep_str):
        from re import compile as recomp
        prog = recomp(grep_str, flags = re.I)
        try:
            while True:
                chan_indx = self.next_chan()
                chan_name = self.chan_name(chan_indx)
                result = prog.search(chan_name)
                if (result):
                    return chan_indx
        except StopIteration:
            return 'name not found'

    def chan_names(self):
        adc_l = map(lambda read: 'adc_' + str(read), np.r_[0:16])
        chans = self.header['multi-chan_inf']['sADCChannelName'][0]
        sampled_chans = self._read_seq()
        #these list of sampled chans is in the order it was sampled
        for num, sampled_chan in enumerate(sampled_chans):
            print '%-3s' '%-3s' '%-8s' '%-10s' '%-10s' %(num, '-'*3,
                                                         adc_l[sampled_chan],
                                                         '-'*8, chans[sampled_chan])

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
        nrows = total_aq/numchans

        # handle optional kwds, for subsetting data
        # start_row and needs to be transulated into byte offset
        # other kwds, num_rows or stop_row do not
        if 'start_row' in kwds.keys():
            start_row = kwds.pop('start_row')
            offset += (start_row * self.base_size)
        if 'num_rows' in kwds.keys():
            nrows = kwds.pop('num_rows')
        if 'stop_row' in kwds.keys():
            # check if start_row is beginning
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

    def get_synch_array(self):
        from abf_header_defs import ABF_BLOCKSIZE
        self.fid.seek(self.header['f_structure']['lSynchArrayPtr'][0] * ABF_BLOCKSIZE)
        synch_array_dtype = [('start', np.int32), ('length', np.int32)]
        synch_array = np.fromfile(self.fid,
                                  synch_array_dtype,
                                  self.header['f_structure']['lSynchArraySize'])
        return synch_array

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

    def start_time(self):
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
            self._File_Time['hour'] = seconds_time/(60*60)
            self._File_Time['minute'] = (seconds_time%(60*60))/60
            self._File_Time['second'] = (seconds_time%(60*60))%60
            self._File_Time['microsecond'] = \
              self.header['environment_inf']['nFileStartMillisecs'][0]\
              * 1000

            #for reading self._File_Time = t_d
            t_d = self._File_Time
            self._file_start_time = datetime(t_d['year'],\
                                      t_d['month'] , t_d['day'],\
                                      t_d['hour'],t_d['minute'],\
                                      t_d['second'],t_d['microsecond'])
            return self._file_start_time
