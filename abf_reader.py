import pdb
import os
import numpy as np
from tempfile import mkdtemp

class waveform(object):
    def __init__(self, **kwds):
        # I think I should make the start index absolute, to the abf file.
        # self._len = 0
        super(waveform, self).__init__(**kwds)

    def __call__(self):
        pass

    def _len(self):
        pass

class waveform_collection(waveform):
    def __init__(self, list_of_waveforms, **kwds):

        # init super first so don't write over self._len
        super(waveform, self).__init__(**kwds)
        self.waveforms = list_of_waveforms

        # make sure each item in list is either wavefrom or
        # waveform_collection type
        for wvf in self.waveforms:
            assert(issubclass(type(wvf),waveform))
        
        # preindex the waveform length
        self._len = 0
        for wvf in self.waveforms:
            self._len += len(wvf())

    def next(self):
        try:
            self._current_waveform += 1
            if (self._current_waveform+1>len(self.waveforms)):
                raise StopIteration            
        except AttributeError:
            self._current_waveform = 0
        return (self.waveforms[self._current_waveform])
    
    def _rewind(self):
        try:
            del self._current_waveform
        except AttributeError:
            pass

    def __call__(self, cmd = '__call__'):
        from numpy import zeros
        length = self._len

        # preindex an array to hold the stitched array
        # call the first wvf to get its shape
        # rewind _current_waveform here
        self._rewind()
        wvf = self.next()
        d = wvf.__getattribute__(cmd)()
        datashp = wvf.__getattribute__(cmd)().shape
        if np.size(datashp)==1:
            out_data = np.zeros((self._len,))
        elif np.size(datashp)==2:
            out_data = np.zeros((self._len, datashp[1]))
        plc_hldr = 0
        self._rewind()

        while True:
            try:
                wvf = self.next()
                td = wvf.__getattribute__(cmd)()
                out_data[plc_hldr:(plc_hldr+len(td))] = td
                plc_hldr += len(td)
            except StopIteration:
                return out_data

# floating methods to be added to abf_waveform on parsing the
# header, e.i. for steps and ramps

def dac_step(self):
    from numpy import tile
    epsd_num = self._epsd_num
    lvl = self._level()
    leng = self._len()
    return tile(lvl, leng)

def dac_ramp(self, epsd_num):
    from numpy import r_
    epsd_num = self._epsd_num
    leng = self._len()
    try:
        return (self._ramp(r_[0:leng]))
    except AttributeError:
        from scipy.interpolate import interp1d
        # make interpolation function, prob better way
        prv = self._prev_epch_lvl()
        lvl = self._level()
        self._ramp = interp1d([0,leng],
                              [prv,lvl])
        return (self._ramp(r_[0:leng]))

class dac_waveform(waveform_collection):
    def __init__(self, abf_reader, dac_num, **kwds):
        assert (dac_num == 0 or dac_num == 1)
        self._dac_num = dac_num
        self._abr = abf_reader
        self._list_o_epochs = []
        
        # construction is complicated involve three methods below,
        # find_actv_epchs, _make_epch, _add_next
        # this could be much cleaner.
        self.find_actv_epchs()
        for i in range(self._len_epch_list):
            self._add_next()
        super(dac_waveform, self).__init__(self._list_o_epochs, **kwds)

    def __call__(self, cmd = '__call__'):
        try:
            return super(dac_waveform, self).__call__(cmd = cmd)
        except AttributeError:
            self.rewind_episode()
            self.next_episode()
            return super(dac_waveform, self).__call__(cmd = cmd)

    def find_actv_epchs(self):
        epch_types = self._abr.header['ext_epoch_waveform_pulses']\
            ['nEpochType'][0, self._dac_num]
        self._actv_epchs = np.nonzero(epch_types)[0]
        self._len_epch_list = len(self._actv_epchs)
        return (self._len_epch_list)

    def _make_epch(self, epch_num):
        eewp = self._abr.header['ext_epoch_waveform_pulses']
        level_init = eewp['fEpochInitLevel'][0]\
            [self._dac_num][epch_num]
        dur_init = eewp['lEpochInitDuration'][0]\
            [self._dac_num][epch_num]
        epoch_type = eewp['nEpochType'][0]\
            [self._dac_num][epch_num]
        level_incrm = eewp['fEpochLevelInc'][0]\
            [self._dac_num][epch_num]
        dur_incrm = eewp['lEpochDurationInc'][0]\
            [self._dac_num][epch_num]
        tmp_epoch = epoch(self._abr, self,
                          epch_num, epoch_type,
                          level_init, dur_init,
                          level_incrm, dur_incrm)
        return tmp_epoch

    def _add_next(self):
        try:
            self._current_epoch += 1
            if (self._current_epoch+1>self._len_epch_list):
                raise StopIteration            
        except AttributeError:
            self._current_epoch = 0

        # make sure to look up the indx in the list of actv epchs,
        # incase their are inactive epochs, this will aviod them
        epch_num = self._actv_epchs[self._current_epoch]
        tmp_epch = self._make_epch(epch_num)
        self.append(tmp_epch)

    def append(self, epoch):
        self._list_o_epochs.append(epoch)

    def get_epoch(self, ordinal):
        assert (ordinal > 0)
        return self._list_o_epochs[ordinal-1]

    def next_episode(self):
        for epch in self._list_o_epochs:
            epch.next()
        return self

    def rewind_episode(self):
        [epch._rewind_episode() for epch in self._list_o_epochs]

class abf_waveform(waveform):
    def __init__(self, abf_reader, start_row, num_rows, **kwds):
        super(abf_waveform, self).__init__(**kwds)
        self._start_row = start_row
        self._len = num_rows
        self._abr = abf_reader

    def __call__(self):
        return self._abr.read_data(start_row = self._start_row,
                                   num_rows = self._len)

class epoch(waveform):
    def __init__(self, abf_reader, dac_wvfrm,
                 ordinal, epoch_type,
                 level_init, dur_init,
                 level_incrm, dur_incrm, **kwds):
        super(epoch, self).__init__(**kwds)
        self._abr = abf_reader
        self._dac_waveform = dac_wvfrm
        self._ordinal = ordinal
        self._type = epoch_type
        self._level_init = level_init
        self._dur_init = dur_init
        self._level_incrm = level_incrm
        self._dur_incrm = dur_incrm
        self.next()

        # dynamically set a level method, based on the epoch type
        import types
        if self._type == 1:
            # type is step
            self.__dict__['level'] = types.MethodType(dac_step, self)
        if self._type == 2:
            # type is ramp
            self.__dict__['level'] = types.MethodType(dac_ramp, self)

        self._th = self._abr.header['trial_hierarchy']
        self._num_episodes = self._th['lEpisodesPerRun'][0]

    def next(self):
        try:
            self._epsd_num += 1
            if (self._epsd_num+1>self._num_episodes):
                raise StopIteration            
        except AttributeError:
            self._epsd_num = 0
        return self

    def __iter__(self):
        return self

    def _set_epsd_num(self, epsd_num):
        self._epsd_num = epsd_num

    def _rewind_episode(self):
        try:
            del self._epsd_num
        except AttributeError:
            pass
        self.next()
            
    def _set_epsd_num(self, epsd_num):
        self._epsd_num = epsd_num

    def _prev_epch_lvl(self):
        if self._ordinal == 0:
            return 0
        else:
            prv_epch = self._dac_waveform.get_epoch(self._ordinal)
            prv_epch._set_epsd_num(self._epsd_num)
            return (prv_epch.level()[-1])

    def get_dp_pad(self):
        SmplsPerEpsd = self._abr.header['trial_hierarchy']\
            ['lNumSamplesPerEpisode'][0]
        RowsPerEpsd = (SmplsPerEpsd/self._abr.num_chans())
        
        # FINALLY FIGURED OUT HOW THE PRE AND POST HOLDS ARE
        # DETERMINED, THEY ARE EPISODE LENGTH / 64 (MODULO DIV)
        # from the pclamp guide on LTP of all things
        dp_one_side_pad = int(RowsPerEpsd) / int(64)
        return dp_one_side_pad

    def _smpls_per_epsd(self):
        SmplsPerEpsd = self._th['lNumSamplesPerEpisode'][0]
        return (SmplsPerEpsd/self._abr.num_chans())

    def _dp_start(self):
        if self._ordinal == 0:
            return ((self._epsd_num * self._smpls_per_epsd()) +
                    self.get_dp_pad())
        else:
            prv_epch = self._dac_waveform.get_epoch(self._ordinal)
            prv_epch._set_epsd_num(self._epsd_num)
            return (prv_epch._dp_end())

    def _dp_end(self):
        if self._ordinal == 0:
            return (self._dp_start() + self._len())
        else:
            prv_epch = self._dac_waveform.get_epoch(self._ordinal)
            prv_epch._set_epsd_num(self._epsd_num)
            return (prv_epch._dp_end() +
                    self._len())

    def _level(self):
        # this returns the _level property, different than the
        # inflated command level
        epsd_num = self._epsd_num
        if self._level_incrm == 0:
            return self._level_init
        else:
            return (self._level_init + (epsd_num * self._level_incrm))

    def _len(self):
        epsd_num = self._epsd_num
        if self._dur_incrm == 0:
            return self._dur_init
        else:
            return (self._dur_init + (epsd_num * self._dur_incrm))

    def __call__(self):
        start_row = self._dp_start()
        num_rows = self._len()
        return (self._abr.read_data(start_row = start_row,
                            num_rows = num_rows))

class episode(abf_waveform):
    def __init__(self, abf_reader, start_row, num_rows, episode_num, **kwds): 
        self._epsd_num = episode_num
        super(episode, self).__init__(abf_reader, start_row, num_rows, **kwds)

class abf_epsd_cnstrctr(object):
    def __init__(self, abf_reader):
        self._abr = abf_reader
        self._th = self._abr.header['trial_hierarchy']
        self._eewp = self._abr.header['ext_epoch_waveform_pulses']
        self._num_episodes = self._th['lEpisodesPerRun'][0]
        
    def _smpls_per_epsd(self):
        SmplsPerEpsd = self._th['lNumSamplesPerEpisode'][0]
        return (SmplsPerEpsd/self._abr.num_chans())

    def __iter__(self):
        return self

    def _num_episodes(self):
        return (self._th['lEpisodesPerRun'])

    def next(self):
        try:
            self._current_epsd += 1
            if (self._current_epsd+1>self._num_episodes):
                raise StopIteration            
        except AttributeError:
            self._current_epsd = 0
        start_row = (self._current_epsd * self._smpls_per_epsd())
        num_rows = self._smpls_per_epsd()
        return episode(self._abr, start_row, num_rows, self._current_epsd+1)

        
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

    def _reopen_editing(self):
        self.fid.close()
        self.fid = file(os.path.join(self.path, self.fname), 'r+b')

    def _reopen_reading(self):
        self.fid.close()
        self.fid = file(os.path.join(self.path, self.fname), 'rb')

    def set_start_time(self, time):
        from math import floor
        self._reopen_editing()
        from datetime import datetime
        # have to create string like yyyymmdd
        assert type(time) == type(datetime.now()), "Only python datetime objects accepted here"
        # have a time date object, good
        yyyymmdd = time.strftime("%Y%m%d")
        lFileStartDate = np.int32(yyyymmdd)

        # 'lFileStartTime is in seconds. find the number of
        # seconds since the start of day
        day_start = datetime(time.year, time.month, time.day)
        FileStartTime = floor((time-day_start).total_seconds())
        lFileStartTime = np.int32(FileStartTime)

        # 'nFileStartMillisecs'
        FileStartMillisecs = ((time-day_start).total_seconds()) % 1 * 1000
        nFileStartMillisecs = np.int16(FileStartMillisecs)

        # change the values of the header, and write entire header back to file
        self.header['fid_size_info']['lFileStartDate'][0] = lFileStartDate
        self.header['fid_size_info']['lFileStartTime'][0] = lFileStartTime
        self.header['environment_inf']['nFileStartMillisecs'][0] = nFileStartMillisecs

        self.fid.seek(0)
        self.fid.write(self.header)
        self.fid.flush()
        self._reopen_reading()

    def stop_watch_time(self):
        return int(self.header['fid_size_info']['lStopwatchTime'][0])
