import pdb
import os
import numpy as np
from tempfile import mkdtemp

class abf_header():
    from abf_header_dtype import abf_header_dtype
    _headr_struct = abf_header_dtype

class abf_reader(object, abf_header):
    def __init__(self, fname):
        import abf_header_defs as defs
        self.defs = defs
        import os
        import numpy as np
        self.fname = os.path.basename(fname)
        self.path = os.path.dirname(os.path.abspath(fname))
        self.path_file = self.path + os.sep + self.fname
        self.fid = file(self.path_file, 'rb')
        self.read_header()
        self.addGain()

    def read_header(self):
        import numpy as np
        self.fid.seek(0)
        self.header = np.fromfile(self.fid, self._headr_struct, 1)

    def _read_seq(self):
        return filter(lambda read: \
            read != -1, self.header['multi-chan_inf']['nADCSamplingSeq'][0])

    def num_chans(self):
        return len(self._read_seq(self))

    def read_data(self):
        '''reads multiplexed data from abfs into an array'''
        ## want to have this method return an iterator that provides data until its out.
        ## have to use data from the header to do this
        ## the times that are asso with discontinuous recording are wonky
        from numpy import fromfile, float32, int16, memmap, float
        offset = self.header['f_structure']['lDataSectionPtr'] * self.defs.ABF_BLOCKSIZE
        total_aq = self.header['fid_size_info']['lActualAcqLength'][0]
        numchans = self.header['trial_hierarchy']['nADCNumChannels'][0]
        howmanyreads = self.header['trial_hierarchy']['lNumSamplesPerEpisode'][0] \
            *  self.header['fid_size_info']['lActualEpisodes'][0] / numchans
        read_seq = self._read_seq()
        ncols = len(read_seq)
        nrows = total_aq/numchans
        reads = map(lambda read: 'adc_' + str(read), read_seq)
        if self.header['f_structure']['nDataFormat'][0]==1: #float data
            data = memmap(self.fid, dtype = float32, shape = (nrows,ncols), offset = offset)
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
        for indx, chan in enumerate(self._read_seq()):
            divis = (self.header['multi-chan_inf']['fInstrumentScaleFactor'][0][chan] * \
                     self.header['multi-chan_inf']['fSignalGain'][0][chan] * \
                     self.header['multi-chan_inf']['fADCProgrammableGain'][0][chan] * \
                     self.addGain[chan])
            mult =  self.header['hardware_inf']['fADCRange'][0] \
                   / self.header['hardware_inf']['lADCResolution'][0]
            offs = self.header['multi-chan_inf']['fInstrumentOffset'][0][chan] - \
                   self.header['multi-chan_inf']['fSignalOffset'][0][chan]
            self.mm[:,indx] = self.mm[:,indx] / divis * mult + offs
        return self.mm
        # inorder to plot this with out decimating data, I will have to
        # 'paginate' this into mutliple files, and then stitch them I should
        # be able to decimate, aswell, especially with data sampled around
        # 10k, with out much affect on appearance also, right now this is

    def addGain(self):
        '''method helps with scaling'''
        self.addGain = self.header['ext_environment_inf']['nTelegraphEnable'][0] * self.header['ext_environment_inf']['fTelegraphAdditGain'][0]
        self.addGain = np.where(self.addGain==0, 1, self.addGain)

    def get_synch_array(self):
        import numpy as np
        self.fid.seek(self.header['f_structure']['lSynchArrayPtr'][0]*self.defs.ABF_BLOCKSIZE)
        synch_array_dtype = [('start', np.int32), ('length', np.int32)]
        synch_array = np.fromfile(self.fid, synch_array_dtype, self.header['f_structure']['lSynchArraySize'])
        return synch_array

    def get_dp_pad(self):
        import numpy as np
        #the total number of data points aquired minus the number of
        #data points defined in the epochs gives the padding data in
        #the file, half of this pad is at the beginning of the file,
        #and half at the end

        #because epochs can be defined as off, meaning no data are
        #aquired, have to filter epochs by epoch type != 0, meaning
        #data are aquired.
        self._actv_wv = np.nonzero(self.header['ext_epoch_waveform_pulses']['nWaveformEnable'][0])[0]
        if len(self._actv_wv)>1:
            raise IndexError('this is currently\
                              not ready to handle\
                              two active waveforms')
        else:
            self._actv_wv = self._actv_wv[0]

        dp_in_epochs = self.header['ext_epoch_waveform_pulses']\
                     ['lEpochInitDuration'][0][self._actv_wv].sum()

                #this is the subtraction
        self._dp_pad =\
        (self.header['trial_hierarchy']['lNumSamplesPerEpisode'][0] \
        / self.header['trial_hierarchy']['nADCNumChannels'][0]) \
        - dp_in_epochs
        #half of this space added before the epoch, and half after, so
        #compute the left pad for convience
        self._dp_lpad = self._dp_pad / 2

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
          return self._sample_rate

    def waveforms_from_header(self):
        import numpy as np

        

    def epochs_from_header(self):
        # have to figure out how to deal with active vs. inactive
        # epochs might be glued together in an object, a 'constructed waveform',

        import numpy as np
        self.episodes_from_header()
        
        #eventually, create an two item list of lists: each sublist
        #will be a list of epoch objects
        #but for now:
        self.epochs = []
        

        #make var for part of header I am pulling the epoch definitions from for
        #readibility
        eewp = self.header['ext_epoch_waveform_pulses']

        #get the number of epochs, by using the type defintions. There are 10
        #possible epochs, those whos type is not zero are actually
        #defined/inuse, not sure if I should use the first [0][0] or the
        #second [0][1] block of the epoch definition, seems to depend on the file.
        self._actv_wv = np.nonzero(eewp['nWaveformEnable'][0])[0]
        if len(self._actv_wv)>1:
            raise IndexError('this is currently not ready to handle to active waveforms')
        else:
            self._actv_wv = self._actv_wv[0]
        self._num_epochs = np.nonzero(eewp['nEpochType'][0][self._actv_wv])[0].size

        # files are longer than protocol for some reason, here is a kludge to
        # get the indexes of each epoch correct, note that I am using the
        # extended epoch wavefrom pulses portion of the header
        try:
            self._dp_pad
        except AttributeError:
            self.get_dp_pad()

        tmp_dp_counter = 0
        tmp_dp_counter += self._dp_lpad

        #iterate over number of epochs
        for epoch in range(self._num_epochs):
            #there is something funny going on with epoch and abf_waveforms, i
            #don't really understand, but to get the data point index for the
            #start of an epoch, I need to use the epoch definitions in the
            #extended epoch epoch waveform pulses section of the header?  this
            #extended header has two of everything, (hence the
            #[0][0][epoch_num] indexing), i guess bc there are two DAC
            #channels that support waveforms?

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
            pdb.set_trace()
            tmp_epoch = abf_epoch(tmp_epoch_level_init,
                          tmp_epoch_dur_init, tmp_epoch_type,
                                  tmp_dp_counter, self,
                                  tmp_epoch_level_incrm,
                                  tmp_epoch_dur_init)
            tmp_dp_counter += tmp_epoch_dur_init
            self.epochs.append(tmp_epoch)

    def episodes_from_header(self):
        '''This creates episodes, but i think episodes composed of epochs with increasing or decreasung durations will break this'''
        import numpy as np
        self.episodes = []
        t_h = self.header['trial_hierarchy']
        self._episode_len = t_h['lNumSamplesPerEpisode'] /\
            t_h['nADCNumChannels']
        self._num_episodes = t_h['lEpisodesPerRun']
        for episode in range(self._num_episodes):
            tmp_episode = abf_episode(self._episode_len * episode, episode)
            self.episodes.append( tmp_episode )
        self.episodes = np.array(self.episodes)

    def lp_filter(self, hertz):
        import mpLOWPASS
        try:
            self.data
            if not len(self.data.shape)==1:
                self.read_data()
        except AttributeError:
            self.read_data()
        for name in self.data.dtype.names:
            self.data[name] = mpLOWPASS.mpLOWPASS(self.data[name], hertz, self.sample_rate())

    def episode_data(self):
        '''reshape continous data into composite episodes, right now incrementing episode/epoch lengths are not supported'''
        self.episodes_from_header()
        self.read_data()
        self.mm = self.mm.reshape((self._num_episodes, self._episode_len, self.mm.shape[1]))
        return self.mm

    def start_time(self):
        try:
            return self._file_start_time
        except AttributeError:
            import datetime
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
            self._file_start_time = datetime.datetime(t_d['year'],\
                                      t_d['month'] , t_d['day'],\
                                      t_d['hour'],t_d['minute'],\
                                      t_d['second'],t_d['microsecond'])
            return self._file_start_time

    def eert_data(self, epoch = 0, episode = 0, run = 0, trail = 0):
        ''' method returns slices of data cooresponding to the designated Epoch Episode Trial Run ### runs and trail are not implemented now'''
        self.episodes_from_header()
        self.epochs_from_header()
        self.episode_data()
        epch_start = self.epochs[epoch]._dp_start
        epch_end = epch_start + self.epochs[epoch]._dur
        return self.mm[episode,epch_start:epch_end,:]

    def chan_names(self):
        import numpy as np
        adc_l = map(lambda read: 'adc_' + str(read), np.r_[0:16])
        chans = self.header['multi-chan_inf']['sADCChannelName'][0]
        sampled_chans = self._read_seq()
        #these list of sampled chans is in the order it was sampled
        for num, sampled_chan in enumerate(sampled_chans):
            print '%-3s' '%-3s' '%-8s' '%-10s' '%-10s' %(num, '-'*3, adc_l[sampled_chan], '-'*8, chans[sampled_chan])

    def chan_name(self, chan_no=0):
        import numpy as np
        chans = self.header['multi-chan_inf']['sADCChannelName'][0]
        sampled_chans = self._read_seq()
        #rstrip removes the trailing white space
        return chans[sampled_chans[chan_no]].rstrip()

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

class abf_waveform:
    # data structure is modified from abf_header_dtype, made singular:
    from abf_header_defs import ABF_EPOCHCOUNT
    dtype = np.dtype([\
        ('nWaveformEnable', np.int16),
        ('nWaveformSource', np.int16),
        ('nInterEpisodeLevel', np.int16),
        ('nEpochType', np.int16, ABF_EPOCHCOUNT),
        ('fEpochInitLevel', np.float32, ABF_EPOCHCOUNT),
        ('fEpochLevelInc', np.float32, ABF_EPOCHCOUNT),
        ('lEpochInitDuration', np.int32, ABF_EPOCHCOUNT),
        ('lEpochDurationInc', np.int32, ABF_EPOCHCOUNT),
        ('nDigitalTrainValue', np.int16, ABF_EPOCHCOUNT), ## 2 * 10 = 20 bytes
        ('nDigitalTrainActiveLogic', np.int16),   ## 2 bytes
        ('sUnused012', np.str_, 18)
        ])

    def __init__(self):
         self._data = np.zeros(1, dtype = self.dtype)
    
    def _from_header(self, abf_reader, which = 0):
         eewp = abf_reader.header['ext_epoch_waveform_pulses']
         for name in eewp.dtype.names:
             # most of these fields are duped for each waveform
             try:
                 self._data[name] = eewp[name][0][which]
             except IndexError:
                 # except nDigitalTrainActiveLogic
                 self._data[name] = eewp[name][0]


    # def __call__(self, episodes_slice = slice(0,9999)):
    #     return_list = []
    #     for episode_indx in range(num_episodes)[episodes_slice]:
    #         for 
            
                     
