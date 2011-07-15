from numpy import memmap, fromfile, float32, int16, memmap, float
from matplotlib.cbook import iterable

import pdb
class abf_chunker(object):
    ''' use this to return data in chunks from abf files'''
    def __init__(self, abr):
        self.abr = abr
        self.dp = self.abr.header['fid_size_info']['lActualAcqLength'][0]
        self._chunk_counter = 0
        self.header_size = self.abr.header['f_structure']['lDataSectionPtr']\
                           * self.abr.defs.ABF_BLOCKSIZE
        self.numchans = self.abr.header['trial_hierarchy']['nADCNumChannels'][0]
        self.howmanyreads = self.abr.header['trial_hierarchy']['lNumSamplesPerEpisode'][0] \
            *  self.abr.header['fid_size_info']['lActualEpisodes'][0] / self.numchans
        self.read_seq = self.abr._read_seq()
        self.ncols = len(self.read_seq)
        self.nrows = self.dp/self.numchans
        self.left = 0
        self.right = self.nrows
        self.set_range((self.left, self.right))
        self.make_chunk_size()

    def set_range(self, left = None, width = None):
        # style points to matplotlib, yay free software
        if width is None and iterable(left):
            if map(type, left)==[int]*2:
                self.left, self.right = left
            elif map(type, left)==[float]*2:
                tested_prcnt_rng = filter(lambda a: a>=0 and a<=1, left)
                if len(tested_prcnt_rng)==2:
                    self.left = int(self.nrows * left[0])
                    self.right = int(self.nrows * left[1])
                else:
                    raise ValueError('if using floats in range, I am \
                    guessing you want to specify and range by percent,\
                    your numbers must be btwn 0, 1')
        if type(width)==int:
            self.left = left
            self.right = left+width
        if type(width)==float and 0<=width<=1:
            self.left = left
            self.right = left + int(self.nrows*width)
        self.range_rows = self.right - self.left

    def make_chunk_size(self, nominal_chunksize = 2**21):
        # make a sane chunksize, so that data can be broken up nicely
        # about 4 megabytes (2**21)? is good?
        if self.abr.header['f_structure']['nDataFormat'][0]==1: #float data
            base_size = 4 * self.ncols # 4byte size per float
        elif self.abr.header['f_structure']['nDataFormat'][0]==0: #integer data
            base_size = 2 * self.ncols # 2byte size per int

        #find number of rows in a 2**21 chunk
        self.chunk_row_size = nominal_chunksize / base_size
        self.num_chunk_rows = self.range_rows / self.chunk_row_size
        self.rmndr_row_size = self.range_rows % self.chunk_row_size
        self.chunksize =  self.chunk_row_size * base_size
        self.num_chunks = self.num_chunk_rows
        self.rmndr_chunk_size = self.rmndr_row_size * base_size
        self.offset_base = self.header_size + (self.left * base_size)
        return

    def __iter__(self):
        return self

    def next(self):
        #check to see if at end of chunks
        if self._chunk_counter==self.num_chunks:
            offset = self._chunk_counter * self.chunksize
            row_size = self.rmndr_row_size
            self._chunk_counter += 1
        elif self._chunk_counter < self.num_chunks:
            offset = self._chunk_counter * self.chunksize
            end_dp = (self._chunk_counter+1) + self.chunksize
            row_size = self.chunk_row_size
            self._chunk_counter += 1
        elif self._chunk_counter > self.num_chunks:
            raise StopIteration

        if self.abr.header['f_structure']['nDataFormat'][0]==1: #float data
            data = memmap(self.abr.fid, dtype = float32, shape = (row_size,self.ncols), offset = offset+self.offset_base)
            return data

        elif self.abr.header['f_structure']['nDataFormat'][0]==0: #integer data
            data = memmap(self.abr.fid, dtype = int16, shape = (row_size,self.ncols),
                          mode = 'r',offset = offset + self.offset_base)
            data = data[:].astype(float32)
            for indx, chan in enumerate(self.read_seq):
                divis = (self.abr.header['multi-chan_inf']['fInstrumentScaleFactor'][0][chan] * \
                         self.abr.header['multi-chan_inf']['fSignalGain'][0][chan] * \
                         self.abr.header['multi-chan_inf']['fADCProgrammableGain'][0][chan] * \
                         self.abr.addGain[chan])
                mult =  self.abr.header['hardware_inf']['fADCRange'][0] / \
                       self.abr.header['hardware_inf']['lADCResolution'][0]
                offs = self.abr.header['multi-chan_inf']['fInstrumentOffset'][0][chan] - \
                       self.abr.header['multi-chan_inf']['fSignalOffset'][0][chan]
                data[:,indx] = data[:,indx] / divis * mult + offs
            return data



