from numpy import memmap, fromfile, float32, int16, memmap, float
from matplotlib.cbook import iterable

import pdb
class abf_chunker(object):
    ''' use this to return data in chunks from abf files'''
    def __init__(self, abr):
        self.abr = abr
        self.dp = self.abr.total_aq
        self._chunk_counter = 0
        self.header_size = self.abr.hdr_offset
        self.numchans = self.abr.num_chans()
        self.ncols = self.numchans
        self.nrows = self.dp/self.numchans

        ###############
        # think like this:
        #       chan1  chan2  chan3
        # row1    1.     0.1    3.
        # row2    0.1    3.     4.

        # make an atomic size, so that data can be broken up with out
        # segmenting cols(channels)
        if self.abr.header['f_structure']['nDataFormat'][0]==1: #float data
            self.base_size = 4 * self.ncols # 4byte size per float
        elif self.abr.header['f_structure']['nDataFormat'][0]==0: #integer data
            self.base_size = 2 * self.ncols # 2byte size per int

        self.left = 0
        self.right = self.nrows
        self.set_range((self.left, self.right))
        self.make_chunk_size()

    def second_to_dp(self, second):
        # second timing ignores num channels of aquisition
        # have to multiply by the row size
        dp = second * self.abr.sample_rate() * self.base_size
        return int(dp)

    def second_to_row(self, second):
        # second timing ignores num channels of aquisition
        # have to multiply by the row size
        row = second * self.abr.sample_rate()
        return int(row)

    def get_range_seconds(self):
        return self.range_rows/self.abr.sample_rate()

    def set_range(self, left = None, width = None, seconds = False):
        '''set the range over the whole file that you want to iterate over, be
        very careful about the types of the range or width argument. 

        If left is a tuple of ints, then these are assumed to be bite ranges

        If left is a tuple of floats btwn 0-1 inclsv, these are assumed to be
        fractions of the x range

        If left is a tuple of floats and seconds is True, these are assumed to
        be seconds of x range'''
        # style points to matplotlib, yay free software

        if width is None and iterable(left):
            if map(type, left)==[int]*2:
                self.left, self.right = left
            elif map(type, left)==[float]*2:
                #if floats are btwn 0-1 inclsv - assume they are fractions of x range
                tested_prcnt_rng = filter(lambda a: a>=0 and a<=1, left)
                if len(tested_prcnt_rng)==2:
                    self.left = int(self.nrows * left[0])
                    self.right = int(self.nrows * left[1])
                #if range is floats and seconds true, 
                elif seconds==True:
                    self.left = self.second_to_row(left[0])
                    self.right = self.second_to_row(left[1])
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

        # set right limit to max, if greater than max.
        if self.right > self.nrows:
            self.right = self.nrows

        # set size of the range of rows STAY IN ROWS (NOT BYTES)
        self.range_rows = (self.right - self.left)
        self.make_chunk_size()

    def make_chunk_size(self, nominal_chunksize = 2**21):

        #local copy of base size, for no reason
        base_size = self.base_size

        #find number of rows in a 2**21 chunk
        self.chunk_row_size = int(nominal_chunksize) / int(base_size)
        self.num_chunk_rows = self.range_rows / self.chunk_row_size
        self.rmndr_row_size = self.range_rows % self.chunk_row_size
        self.chunksize =  self.chunk_row_size * base_size
        self.num_chunks = self.num_chunk_rows
        self.rmndr_chunk_size = self.rmndr_row_size * base_size
        self.offset_base = self.header_size + (self.left * base_size)

        #find the number of chunks, as a float!
        self.frct_chunks = self.num_chunk_rows + (float(self.rmndr_row_size) / self.chunk_row_size)
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
            try:
                data = memmap(self.abr.fid, dtype = int16, shape = (row_size,self.ncols),
                              mode = 'r',offset = offset + self.offset_base)
            except ValueError:
                pdb.set_trace()
            data = data[:].astype(float32)
            data = self.abr.scale_int_data(data)
            return data



