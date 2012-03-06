import os

def rename_abf(fname):
        fname = os.path.basename(fname)
        path = os.path.dirname(os.path.abspath(fname))
        path_file = os.path.join(path, fname)
        fnoext = fname.split(os.extsep)[0]
        ## pull out the year, month and day from the fname
        year = fnoext[0:2]
        month = fnoext[2]
        day = fnoext[3:5]
        fnum = fnoext[5:8]
        ## make yyyy_mm_dd_NNNN
        newyear = '20' + year

        ## because the month is 1 : 9 for jan - sep, than o = october,
        ## n = november, d = decemeber, have to have some cases here


        try:
            int(month)
            newmonth = '0' + month
        except ValueError:
            if month=='o':
                newmonth = '10'
            elif month=='n':
                newmonth = '11'
            elif month == 'd':
                newmonth = '12'
        

        ## day stays same
        ## fnum needs a prepended 0
                
        newfnum = '0' + fnum

        new_name = newyear+'_'+newmonth+'_'+day+'_'+newfnum+os.extsep+'abf'
        new_path_file = os.path.join(path, new_name)
        os.rename(path_file, new_path_file)
