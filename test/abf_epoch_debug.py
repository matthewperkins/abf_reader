# want to read the header into a structure, right from file, not a dictionary.
import pdb

class abf_header():
    from abf_header_dtype import abf_header_dtype
    _headr_struct = abf_header_dtype
        
# thinking about adding a base class, for abf data - this class would always assosiate a file name with what ever data it contains, and ideally a cell name too

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

    def read_header(self):
        import numpy as np
        self.fid.seek(0)
        self.header = np.fromfile(self.fid, self._headr_struct, 1)
        
    def read_data(self):
        '''reads multiplexed data from abfs into an array'''
        ## have to use data from the header to do this
        ## the times that are asso with discontinuous recording are wonky
        import numpy as np
        self.fid.seek(self.header['f_structure']['lDataSectionPtr'] * self.defs.ABF_BLOCKSIZE)
        total_aq = self.header['fid_size_info']['lActualAcqLength'][0]
        numchans = self.header['trial_hierarchy']['nADCNumChannels'][0]
        howmanyreads = self.header['trial_hierarchy']['lNumSamplesPerEpisode'][0] \
            *  self.header['fid_size_info']['lActualEpisodes'][0] / numchans
        read_seq = filter(lambda read: \
            read != -1, self.header['multi-chan_inf']['nADCSamplingSeq'][0])
        reads = map(lambda read: 'adc_' + str(read), read_seq)
        if self.header['f_structure']['nDataFormat'][0]==1: #float data
            data_dtype = zip(reads, [np.float32]*len(reads))
            scaled_data = np.fromfile(self.fid, data_dtype, howmanyreads)
        elif self.header['f_structure']['nDataFormat'][0]==0: #integer data
            data_dtype = zip(reads, [np.int16]*len(reads))
            scaled_data_dtype = zip(reads, [np.float32]*len(reads))
            data = np.fromfile(self.fid, data_dtype, howmanyreads)
            scaled_data = data.astype(scaled_data_dtype)
            for read in reads:
                scaled_data[read] = 1000 * scaled_data[read]/(np.float(self.header['hardware_inf']['lADCResolution'][0]))
        time = np.r_[0:total_aq]
        synch_array = self.get_synch_array()
        for i, row in enumerate(synch_array):
            time[i*row[1]:i*row[1]+row[1]] = time[i*row[1]:i*row[1]+row[1]] + row[0]
        time = time*self.header['trial_hierarchy']['fSynchTimeUnit'][0]/1000000.
        time.dtype = zip(reads, [time.dtype]*len(reads))
        self.data = scaled_data
        self.data_time = time
        return scaled_data

    def get_synch_array(self):
        import numpy as np
        self.fid.seek(self.header['f_structure']['lSynchArrayPtr'][0]*self.defs.ABF_BLOCKSIZE)
        synch_array_dtype = [('start', np.int32), ('length', np.int32)]
        synch_array = np.fromfile(self.fid, synch_array_dtype, self.header['f_structure']['lSynchArraySize'])
        return synch_array

    def get_dp_pad(self):
        self._dp_pad =\
        (self.header['trial_hierarchy']['lNumSamplesPerEpisode'][0] \
        / self.header['trial_hierarchy']['nADCNumChannels'][0]) \
        - self.header['ext_epoch_waveform_pulses']\
                     ['lEpochInitDuration'][0].sum()

        #half of this space added before the epoch, and half after.
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

    def epochs_from_header(self):
        import numpy as np

        #create an empty list to hold epoch objects
        self.epochs = []
        
        #make var for part of header I am pulling the epoch definitions from for
        #readibility
        eewp = self.header['ext_epoch_waveform_pulses']

        #get the number of epochs, by using the type defintions. There are 10
        #possible epochs, those whos type is not zero are actually
        #defined/inuse
        self._num_epochs = np.nonzero(eewp['nEpochType'][0][0])[0].size

        
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

            tmp_epoch_level_init = eewp['fEpochInitLevel'][0][0][epoch]
            tmp_epoch_dur_init = eewp['lEpochInitDuration'][0][0][epoch]
            tmp_epoch_type = eewp['nEpochType'][0][0][epoch]
            tmp_epoch_level_incrm = eewp['fEpochLevelInc'][0][0][epoch]
            tmp_epoch_dur_incrm = eewp['lEpochDurationInc'][0][0][epoch]
            tmp_epoch = abf_epoch(tmp_epoch_level_init, \
                          tmp_epoch_dur_init, tmp_epoch_type, tmp_dp_counter,\
                          tmp_epoch_level_incrm, tmp_epoch_dur_init)
            tmp_dp_counter += tmp_epoch_dur_init
            self.epochs.append(tmp_epoch)

    def episodes_from_header(self):
        self.episodes = []
        t_h = self.header['trial_hierarchy']
        self._episode_len = t_h['lNumSamplesPerEpisode'] /\
            t_h['nADCNumChannels']
        self._num_episodes = t_h['lEpisodesPerRun']
        for episode in range(self._num_episodes):
            tmp_episode = abf_episode(self._episode_len * episode, episode)
            self.episodes.append( tmp_episode )

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
        try:
            self.data
        except AttributeError:
            self.read_data()
        self.data = self.data.reshape(self._num_episodes, self._episode_len)

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
        self.episode_data()
        epch_start = self.epochs[epoch]._dp_start
        epch_end = epch_start + self.epochs[epoch]._dur
        return self.data[episode,epch_start:epch_end]

        ### Fundamental unit is an epidsode, many episodes make a run,
        ### many runs make a trial.

class abf_episode:

    def __init__(self, data_point_start, episode_num):
        # this data point will be relative to the run (runs are not
        # implemented yet, because, all my files / protocols only have a
        # single run.
        self._dp_start = data_point_start
        self._episode_num = episode_num

class abf_epoch:

    def __init__(self, level_init,  dur_init, epoch_type, data_point_start,\
                     level_incrm=0, dur_incrm=0):
        # data point_start is relative to the start of the episode, it is not
        # an index for data
        ### incrementing durations are not supported now.
        self._level_init = level_init
        self._dur_init = dur_init
        self._level_incrm = level_incrm
        self._dur_incrm = dur_incrm
        self._type = epoch_type
        self._level = level_init
        self._dur = dur_init
        self._dp_start = data_point_start
        if epoch_type==1:
            self._type = 'step'
        elif epoch_type==2:
            self._type = 'ramp'

    def __call__(self, abf_episode):
        self._level = self._level_init + (self._level_incrm * abf_episode.episode_num)
        # self._dur = self._dur_init + (self._dur_incrm * iter_num)
        
        # I have not tried any of the incrementing duration protocols, and I
        # am not sure if the episode length will stay the same for all
        # episode, despite the increasing duration of an epoch, for now, I
        # will assume that the episode length in the header is fixed, and so
        # the epochs class should have start and stop data point indexes as
        # properties / fields that are relative to the episode.
        ########## for now, just right this for non- incrementing durations -
        ########## add later if needed

import os
fname = os.path.expanduser('~') + os.sep + '0201-0191.abf'
fname2 = os.path.expanduser('~') + '/Downloads/testabfs/2011_05_11_0069.abf'
#abr = abf_reader(fname)
abr2 = abf_reader(fname2)
#abr.epochs_from_header()
abr2.epochs_from_header()
abr2.episodes_from_header()
