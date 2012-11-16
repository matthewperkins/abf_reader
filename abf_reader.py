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

def parse_ep_stim(abf, DAC_num):
    abf.read_header()
    abf_header = abf.header

    # values to compute
    # sweep_indxs
    # episode_indxs
    # epoch_indxs
    # epoch_levels
    # epoch_types
    
    num_epsds = abf_header['trial_hierarchy']['lEpisodesPerRun'][0]
    actv_epchs = find_actv_epchs(abf_header, DAC_num)
    num_epchs = len(actv_epchs)
    DAC_len = abf.total_aq() / abf.num_chans()
    sweep_len = DAC_len/num_epsds
    pad = get_dp_pad(abf)
    epsd_len = sweep_len - 2*pad

    eewp = abf_header['ext_epoch_waveform_pulses']

    # epoch levels
    level_inits = eewp['fEpochInitLevel'][0]\
        [DAC_num][actv_epchs]
    level_incrms = eewp['fEpochLevelInc'][0]\
        [DAC_num][actv_epchs]
    range_array = np.transpose(np.tile(np.c_[0:num_epsds], (num_epchs,1)))
    tmp_lvl_incrms = np.tile(level_incrms, (num_epsds,1)) * range_array
    levels = np.tile(level_inits,(num_epsds,1))+tmp_lvl_incrms

    # epoch durs
    dur_inits = eewp['lEpochInitDuration'][0]\
        [DAC_num][actv_epchs]
    dur_incrms = eewp['lEpochDurationInc'][0]\
        [DAC_num][actv_epchs]
    tmp_dur_incrms = np.tile(dur_incrms, (num_epsds,1)) * range_array
    durs = np.tile(dur_inits,(num_epsds,1))+tmp_dur_incrms

    # epoch types
    epoch_types = eewp['nEpochType'][0]\
        [DAC_num][actv_epchs]

    # make sweeps slices
    lft = np.r_[0:DAC_len+1:sweep_len][:-1]
    rght = np.r_[0:DAC_len+1:sweep_len][1:]
    sweep_indxs = c_[lft,rght]
    sweep_slc_array = np.array([slice(l,r) for l,r in zip(lft, rght)])

    # make epsd slices
    lft += pad
    rght = lft+epsd_len
    epsd_indxs = c_[lft,rght]
    epsd_slc_array = np.array([slice(l,r) for l,r in zip(lft, rght)])

    #### epoch starts (rel to episode start) ####
    
    # assume left points are fixed in
    # duration changing epochs. I hope this means the maximum duration
    # for an epoch (over the whole trail) must pass before next epoch
    # statrs
    max_dur = np.max(durs,0)
    max_durs = np.tile(max_dur, (num_epsds,1))

    # accumulate durations (add a zero column) to get starts
    strts = c_[repeat(0,num_epsds), max_durs[:,0:-1]]
    epch_strt_indxs = np.add.accumulate(strts,1)

    #### epoch starts (rel to 0) ####
    for i_epch in range(num_epchs):
        epch_strt_indxs[:,i_epch]+=epsd_indxs[:,0]

    epch_end_indxs = epch_strt_indxs+durs

    epch_indx_list = []
    epch_slc_list = []
    for i_epsd in range(num_epsds):
        tmp_epsd_list=[]
        tmp_slc_list=[]
        for i_epch in range(num_epchs):
            t = [epch_strt_indxs[i_epsd,i_epch],
                 epch_end_indxs[i_epsd,i_epch]]
            tmp_epsd_list.append(t)
            tmp_slc_list.append(slice(t[0],t[1]))
        epch_indx_list.append(tmp_epsd_list)
        epch_slc_list.append(tmp_slc_list)
    epch_indxs = np.array(epch_indx_list)

    return (sweep_indxs, epsd_indxs, epch_indxs)

def get_dp_pad(abf):
    SmplsPerEpsd = abf.header['trial_hierarchy']\
        ['lNumSamplesPerEpisode'][0]
    RowsPerEpsd = (SmplsPerEpsd/abf.num_chans())

    # FINALLY FIGURED OUT HOW THE PRE AND POST HOLDS ARE
    # DETERMINED, THEY ARE EPISODE LENGTH / 64 (MODULO DIV)
    # from the pclamp guide on LTP of all things
    dp_one_side_pad = int(RowsPerEpsd) / int(64)
    return dp_one_side_pad

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
        self.addGain = self.header['ext_environment_inf']['nTelegraphEnable'][0] * self.header['ext_environment_inf']['fTelegraphAdditGain'][0]
        self.addGain = np.where(self.addGain==0, 1, self.addGain)

    def get_synch_array(self):
        from abf_header_defs import ABF_BLOCKSIZE
        self.fid.seek(self.header['f_structure']['lSynchArrayPtr'][0] * ABF_BLOCKSIZE)
        synch_array_dtype = [('start', np.int32), ('length', np.int32)]
        synch_array = np.fromfile(self.fid, synch_array_dtype, self.header['f_structure']['lSynchArraySize'])
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

