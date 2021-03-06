from numpy import memmap, fromfile, float32, int16, memmap, issubdtype, all, array
from matplotlib.cbook import iterable

import pdb
class abf_chunker(object):
    ''' use this to return data in chunks from abf files'''
    def __init__(self, abr, **kwds):
        self.abr = abr
        self.dp = self.abr.total_aq()
        self._chunk_counter = 0
        self.header_size = self.abr.hdr_offset()
        self.numchans = self.abr.num_chans()
        self.ncols = self.numchans
        self.nrows = self.dp//self.numchans
        if 'nominal_chunksize' in kwds:
            self._nominal_chunksize = kwds.pop('nominal_chunksize')
        else:
            self._nominal_chunksize = 2**20

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
        super(abf_chunker, self).__init__(**kwds)

    def rewind(self):
        self._chunk_counter = 0

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
        be seconds of x range
        '''
        # style points to matplotlib, yay free software
        if width is None and iterable(left):
            # this needs some clean up
            if all([issubdtype(type(l), int) for l in left]):
                self.left, self.right = left
            elif all([issubdtype(type(l),float) for l in left]):
                #if floats are btwn 0-1 inclsv - assume they are
                #fractions of x range
                al = array(left)
                if seconds==True:
                    # check to make sure the bounds in seconds are pos.
                    # neg values will read the whol file
                    # assert self.abr.get_synch_array().size==0, "File has been paused, not implemented" 
                    for sec in left: assert sec>=0, "seconds bounds must be greater than zero"
                    self.left = self.second_to_row(left[0])
                    self.right = self.second_to_row(left[1])
                elif all((al>=0) & (al<=1)):
                    self.left = int(self.nrows * left[0])
                    self.right = int(self.nrows * left[1])
                #if range is floats and seconds true, 
                else:
                    raise ValueError('if using floats in range, I am \
                    guessing you want to specify and range by percent,\
                    your numbers must be btwn 0, 1')
        elif type(width)==int and type(left)==int:
            assert seconds!=True, 'getting mixed messages! are these ranges seconds? or what?'
            self.left = left
            self.right = left+width
        elif type(width)==float and type(left)==float:
            assert seconds==True, 'getting mixed messages! are these ranges seconds or what?'
            self.left = self.second_to_row(left)
            self.right = self.second_to_row(left+width)
        else:
            raise TypeError('Just use floats you asswhole')

        # set right limit to max, if greater than max.
        if self.right > self.nrows:
            self.right = self.nrows

        # set size of the range of rows STAY IN ROWS (NOT BYTES)
        self.range_rows = int(self.right - self.left)
        self.make_chunk_size()

    def make_chunk_size(self):

        #be default this will make nominal chunks of size
        # 2**27, this can be changed by setting the _nominal_chunksize property

        #local copy of base size, for no reason
        base_size = self.base_size

        #find number of rows in a 2**21 chunk
        self.chunk_row_size = self._nominal_chunksize // base_size
        #if self.range_rows < self.chunk_row_size:
            
        self.num_chunk_rows = self.range_rows // self.chunk_row_size
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

    def __next__(self):
        #check to see if at end of chunks
        if self._chunk_counter==self.num_chunks:
            offset = int(self._chunk_counter * self.chunksize)
            row_size = self.rmndr_row_size
            self._chunk_counter += 1
        elif self._chunk_counter < self.num_chunks:
            offset = int(self._chunk_counter * self.chunksize)
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



