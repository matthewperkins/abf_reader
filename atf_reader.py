class ATFreader(object):
    def __init__(self, path):
        self.fid = file(path, 'rt')
        self._parse_header()
        super(ATFreader,self).__init__()

    def _parse_header(self):
        self.fid.seek(0)
        (ATF,v) = self.fid.readline().rstrip().split('\t')
        assert ATF=='ATF', "not right file type"
        assert v=='1.0', "not right version"
        (NOptRcrds,NColumns) = map(int, self.fid.readline().rstrip().split('\t'))
        self.NOptRcrds = NOptRcrds
        self.NColumns = NColumns

    def read_data(self):
        from numpy import loadtxt
        self.fid.seek(0)
        return loadtxt(self.fid, delimiter='\t', skiprows = self.NColumns+1)
