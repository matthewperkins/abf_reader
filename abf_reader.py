import pdb
import os
import numpy as np
from tempfile import mkdtemp

# I have to decide if I want to use relative or absolute indxs I think
# I should use absolute?  waveform should be generic, having the
# signatures of a generic data source.

# I think I am going to make this class only for deriving
class waveform(object):
    def __init__(self, **kwds):
        # I think I should make the start index absolute, to the abf file.
        self._start_indx = 0
        self._end_indx = 0
        self._len = 0
        super(waveform, self).__init__(**kwds)

    def __call__(self):
        pass

class epoch(waveform):
    def __init__(self, data, **kwds):
        # I think I should make the start index absolute, to the abf file.
        self._data = data
        super(epoch, self).__init__(**kwds)
        self._start_indx = 0
        self._end_indx = len(data)
        self._len = len(data)

    def __call__(self):
        return self._data

class waveform_collection(waveform):
    def __init__(self, list_of_waveforms, **kwds):

        # init super first so don't write over self._len
        super(waveform, self).__init__(**kwds)
        self.waveforms = list_of_waveforms

        # make sure each item is list is either wavefrom or
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
        self._current_waveform = -1

    def __call__(self):
        from numpy import zeros
        out_data = np.zeros(self._len)
        plc_hldr = 0

        # have to rewind _current_waveform here
        self._rewind()
        while True:
            try:
                wvf = self.next()
                td = wvf()
                out_data[plc_hldr:(plc_hldr+len(td))] = td
                plc_hldr += len(td)
            except StopIteration:
                return out_data
    
class eert(object):
    def __init__(self, abf_reader, **kwds):
        self._abr = abf_reader
        parse_header(self)
        super(eert, self).__init__(**kwds)

    def parse_header(self):
        return

    def get_dp_pad(self, waveform = 0):
        '''get dp_pad for the designated  waveform'''
        #the total number of data points aquired minus the number of
        #data points defined in the epochs gives the padding data in
        #the file, half of this pad is at the beginning of the file,
        #and half at the end

        #because epochs can be defined as off, meaning no data are
        #aquired, have to filter epochs by epoch type != 0, meaning
        #data are aquired.
        if ( (waveform < 0) | (waveform > 1)):
            raise IndexError('waveform must be 0 or 1')
        self._actv_wv = \
            np.nonzero(self._abr.header['ext_epoch_waveform_pulses']\
                           ['nWaveformEnable'][0, waveform])[0]
        dp_in_epochs = self._abr.header['ext_epoch_waveform_pulses']\
                     ['lEpochInitDuration'][0][self._actv_wv].sum()

        #this is the subtraction
        self._dp_pad =\
        (self._abr.header['trial_hierarchy']['lNumSamplesPerEpisode'][0] \
        / self._abr.header['trial_hierarchy']['nADCNumChannels'][0]) \
        - dp_in_epochs
        #half of this space added before the epoch, and half after, so
        #compute the left pad for convience
        self._dp_lpad = self._dp_pad / 2

    def epochs_from_header(self):
        # have to figure out how to deal with active vs. inactive
        # epochs might be glued together in an object, a 'constructed waveform',

        self.episodes_from_header()
        
        #eventually, create an two item list of lists: each sublist
        #will be a list of epoch objects
        #but for now:
        self.epochs = []
        
        #for readibility
        eewp = self.header['ext_epoch_waveform_pulses']

        self._actv_wv = np.nonzero(eewp['nWaveformEnable'][0])[0]
        if len(self._actv_wv)>1:
            raise IndexError('this is currently not ready to handle to active waveforms')
        else:
            self._actv_wv = self._actv_wv[0]

        #get the number of epochs, by using the type. There are 10
        #possible epochs, those whos type is not zero are actually
        #defined/inuse
        self._num_epochs = np.nonzero(eewp['nEpochType'][0][self._actv_wv])[0].size

        # more data are recorded than the num data points defined in
        # the epochs.

        # Here is a kludge to get the indexes of each
        # epoch correct, note that I am using the extended epoch
        # wavefrom pulses portion of the header
        try:
            self._dp_pad
        except AttributeError:
            self.get_dp_pad()

        tmp_dp_counter = 0
        tmp_dp_counter += self._dp_lpad

        # Make an epoch object for each epoch.
        # At the moment only support abfs with a single active waveform
        for epoch in range(self._num_epochs):
            tmp_epoch_level_init = eewp['fEpochInitLevel'][0]\
                [self._actv_wv][epoch]
            tmp_epoch_dur_init = eewp['lEpochInitDuration'][0]\
                [self._actv_wv][epoch]
            tmp_epoch_type = eewp['nEpochType'][0]\
                [self._actv_wv][epoch]
            tmp_epoch_level_incrm = eewp['fEpochLevelInc'][0]\
                [self._actv_wv][epoch]
            tmp_epoch_dur_incrm = eewp['lEpochDurationInc'][0]\
                [self._actv_wv][epoch]
            tmp_epoch = abf_epoch(tmp_epoch_level_init,
                          tmp_epoch_dur_init, tmp_epoch_type,
                                  tmp_dp_counter, self,
                                  tmp_epoch_level_incrm,
                                  tmp_epoch_dur_init)
            tmp_dp_counter += tmp_epoch_dur_init
            self.epochs.append(tmp_epoch)

    def episodes_from_header(self):
        '''This creates episodes, but i think episodes composed of epochs with increasing or decreasung durations will break this'''
        self.episodes = []
        t_h = self.header['trial_hierarchy']
        self._episode_len = t_h['lNumSamplesPerEpisode'] /\
            t_h['nADCNumChannels']
        self._num_episodes = t_h['lEpisodesPerRun']
        for episode in range(self._num_episodes):
            tmp_episode = abf_episode(self._episode_len * episode, episode)
            self.episodes.append( tmp_episode )
        self.episodes = np.array(self.episodes)

    def episode_data(self):
        '''reshape continous data into composite episodes, right now incrementing episode/epoch lengths are not supported'''
        self.episodes_from_header()
        self.read_data()
        self.mm = self.mm.reshape((self._num_episodes, self._episode_len, self.mm.shape[1]))
        return self.mm

class abf_reader(object):
    def __init__(self, fname):
        from abf_header_dtype import abf_header_dtype
        self._headr_struct = abf_header_dtype
        self.fname = os.path.basename(fname)
        self.path = os.path.dirname(os.path.abspath(fname))
        self.path_file = self.path + os.sep + self.fname
        self.fid = file(self.path_file, 'rb')
        self.read_header()
        self.addGain()
        self._chan_holder = -1

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

    def read_data(self):
        '''reads multiplexed data from abfs into an array'''
        ## want to have this method return an iterator that provides data until its out.
        ## have to use data from the header to do this
        ## the times that are asso with discontinuous recording are wonky
        from numpy import fromfile, float32, int16, memmap, float
        from abf_header_defs import ABF_BLOCKSIZE
        offset = self.header['f_structure']['lDataSectionPtr'] * ABF_BLOCKSIZE
        total_aq = self.header['fid_size_info']['lActualAcqLength'][0]
        numchans = self.header['trial_hierarchy']['nADCNumChannels'][0]
        howmanyreads = self.header['trial_hierarchy']['lNumSamplesPerEpisode'][0] \
            *  self.header['fid_size_info']['lActualEpisodes'][0] / numchans
        read_seq = self._read_seq()
        ncols = len(read_seq)
        nrows = total_aq/numchans
        reads = map(lambda read: 'adc_' + str(read), read_seq)

        #see if is float data
        if self.header['f_structure']['nDataFormat'][0]==1: 
            data = memmap(self.fid, dtype = float32, shape = (nrows,ncols), offset = offset)

        #see if is integer data
        elif self.header['f_structure']['nDataFormat'][0]==0: #integer data
            data = memmap(self.fid, dtype = int16, shape = (nrows,ncols), mode = 'r',offset = offset)

        # make a writeable temporarory memory map for helping? with memusage?,
        tmp_dir = mkdtemp()
        tmp_map_name = os.path.join(tmp_dir, '_tmp_' + self.fname.split(os.extsep)[0])
        self.mm = memmap(tmp_map_name, dtype = float32, shape = (nrows,ncols), mode = 'w+')
        self.mm[:] = data[:].astype(np.float32)

        # delete the memory map to flush from file
        del self.mm

        # then reload?
        self.mm = memmap(tmp_map_name, dtype = float32, shape = (nrows,ncols), mode = 'r+')

        # now do some scaling
        self.mm = self.scale_int_data(self.mm)
        return self.mm

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
          return self._sample_rate[0]

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

class abf_episode:

    def __init__(self, data_point_start, episode_num):
        # this data point will be relative to the run (runs are not
        # implemented yet, because, all my files / protocols only have a
        # single run.
        self._dp_start = data_point_start
        self._episode_num = episode_num

class abf_epoch:

    def __init__(self, level_init,  dur_init, epoch_type, data_point_start,\
                     abf_reader, level_incrm=0, dur_incrm=0):
        # data point_start is relative to the start of the episode, it is not
        # an index for data
        ### incrementing durations are not supported now.
        self._abf_reader = abf_reader
        self._level_init = level_init
        self._dur_init = dur_init
        self._level_incrm = level_incrm
        self._dur_incrm = dur_incrm
        self._type = epoch_type
        self._level = level_init
        self._levels = self.levels()
        self._dur = dur_init
        self._dp_start = data_point_start
        self._dp_end = self._dp_start + self._dur
        if epoch_type==1:
            self._type = 'step'
        elif epoch_type==2:
            self._type = 'ramp'

    def __call__(self):
        levels = self._levels
        return (levels, self.data())

    def levels(self):
        ###NOTE THIS METHOD IS BROKEN, NEEDS TO LOOK FIRST AT THE
        ###['ext_epoch_waveform_pulses']['nWaveformEnable'] and only
        ###return the waveform pulses that are enabled. - ei the DAC
        ###that are out putting stuff in the protocol
        self._levels = []
        for epi_num, episode in enumerate(self._abf_reader.episodes):
            self._levels.append(self._level_init + (self._level_incrm * epi_num))
        self._levels = np.array(self._levels)
        return self._levels

    def data(self):
        self._abf_reader.episode_data()
        # indexs are [episode,datapoints,chan_no]
        return self._abf_reader.mm[:, self._dp_start : self._dp_start + self._dur_init, :]

#    def mk_durs(self):
        # self._dur = self._dur_init + (self._dur_incrm * iter_num)
        # I have not tried any of the incrementing duration protocols, and I
        # am not sure if the episode length will stay the same for all
        # episode, despite the increasing duration of an epoch, for now, I
        # will assume that the episode length in the header is fixed, and so
        # the epochs class should have start and stop data point indexes as
        # properties / fields that are relative to the episode.
        ########## for now, just right this for non- incrementing durations -
        ########## add later if needed

                     
