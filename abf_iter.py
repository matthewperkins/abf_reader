from abf_reader import *
import matplotlib.pyplot as plt

class epsd_iter(object):
    def __init__(self, abf, **kwds):
        self.abf = abf
        self.epsd_indxs = make_epsd_indxs(self.abf.header)
        self._num_epsds = get_num_episodes(self.abf.header)
        super(epsd_iter, self).__init__(**kwds)

    def __iter__(self):
        return self.gen_iter()

    def gen_iter(self):
        # using generator
        for index in range(self._num_epsds):
            yield self.abf.read_data(\
                start_row = self.epsd_indxs[index,0],
                stop_row = self.epsd_indxs[index,1])

if __name__=='__main__':
    ''' debug script '''
    import os
    labdir = os.environ.get("LABDIR")
    abfpath = os.path.join(labdir,
                 'B44_B48_experiments',
                 'b48_fI_summary',
                 'frf', '2012_08_24_0004.abf')
    abf = abf_reader(abfpath)
    test_epsd_iter = epsd_iter(abf)
    for epsd in test_epsd_iter:
        plt.plot(epsd[:,0], '-g', linewidth = 0.4)
    plt.show()
    
        

        
        

    

    
