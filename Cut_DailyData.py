import pickle
from obspy.core import Stream, AttribDict
from obspy import UTCDateTime, read
from os.path import join,basename
import pandas as pd
import sys
import obspy
import subprocess
import os


def process(db, network, station):

    sta_time = UTCDateTime("2018-01-01T00:00:00")
    end_time = UTCDateTime("2019-01-01T00:00:00")
    data_mini_length = 1e4  ## Minimal data length (seconds)


    key = '{}.{}'.format(network, station)
    sta = db[key]

    old_dir = 'Raw_Data'
    new_dir = 'Selected_Data'

    ts = sta_time
    te = ts + 3600 * 24

    while ts < end_time:
    
        y_str = str(ts.year)
        jul_str = str(ts.julday).zfill(3)

        file_name = ['SAC_PZ',sta.network,sta.station,y_str]
        file_name = '.'.join(file_name)
        file_path = [sta.station,old_dir,file_name]
        file_path = '/'.join(file_path)

        try:

            name_d = [file_path,jul_str,'BHZ.SAC']
            name_d = '.'.join(name_d)
            name_p = [file_path,jul_str,'BPR.SAC']
            name_p = '.'.join(name_p)
            sac_d = read(pathname_or_url = name_d, format="SAC")
            sac_p = read(pathname_or_url = name_p, format="SAC")

        except Exception:

            ts += 3600 * 24
            te = ts + 3600 * 24
            print("Failed on " + y_str + '.' + jul_str)
            continue
    

        filename_prefix = [sta.station,new_dir,y_str + '.' + jul_str]
        filename_prefix = '/'.join(filename_prefix)
        print("working on " + filename_prefix)

        disp = sac_d.slice(ts, te)
        pres = sac_p.slice(ts, te)
        print('Done slice')

        if len(disp.traces) == 0 or len(pres.traces) == 0:

            ts += 3600 * 24
            te = ts + 3600 * 24
            continue

        else:

            # Get the pressure and displacement are the same length
            try:

                disp, pres = unify_data_length(disp, pres, data_mini_length)
                print("Done unify the length")

            except AttributeError:

                print("Caught")
                ts += 3600 * 24
                te = ts + 3600 * 24
                continue

        fileZ = filename_prefix + ".BHZ"
        fileP = filename_prefix + ".P"

        trZ = disp.traces[0]
        trP = pres.traces[0]

        trZ.write(fileZ, format='sac')
        trP.write(fileP, format='sac')

        time_correction(fileZ, fileP)

        ts += 3600 * 24
        te = ts + 3600 * 24


def time_correction(name_BHZ, name_BPR):

    stream_BPR = obspy.read(name_BPR)
    trace_BPR = stream_BPR[0]
    starttime_BPR = trace_BPR.stats.starttime
    endtime_BPR = trace_BPR.stats.endtime

    stream_BHZ = obspy.read(name_BHZ)
    trace_BHZ = stream_BHZ[0]
    starttime_BHZ = trace_BHZ.stats.starttime
    endtime_BHZ = trace_BHZ.stats.endtime

    starttime = max(starttime_BPR, starttime_BHZ)
    endtime = min(endtime_BPR, endtime_BHZ)

    daytime = obspy.UTCDateTime(starttime.year, starttime.month, starttime.day, 0, 0, 0, 0)
    btime = starttime - daytime
    etime = endtime - daytime

    if starttime < endtime:

        trace_BHZ.trim(starttime, endtime)
        trace_BPR.trim(starttime, endtime)
        stream_BHZ.write(name_BHZ, format='sac')
        stream_BPR.write(name_BPR, format='sac')

        p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        s = ''
        s += 'r %s %s\n' % (name_BHZ, name_BPR)
        s += 'ch nzyear %d nzjday %d nzhour %d nzmin %d nzsec %d nzmsec %d\n' % (starttime.year, starttime.julday, 0, 0, 0, 0)
        s += 'ch b %.2f e %.2f\n' % (btime, etime)
        s += 'wh\n'
        s += 'q\n' 
        p.communicate(s.encode())


def unify_data_length(d, p, mini_length):

    """
    Find the time and return the new Stream

    :param d,p: a Stream object has several time slice
    :return d_slice,p_slice: Stream object (For remove response)
    :exception pressure date not coherent with displacement
                data length too short

    """

    d = d.traces
    index = 0
    tempT = d[0]

    for i in range(len(d)):
    
        if len(d[i]) > len(tempT):
            index = i
            tempT = d[i]

    d_slice = d[index]

    #  the admittance of each day is based on 5 * 2000 s long time series

    if d_slice.meta.endtime - d_slice.meta.starttime < mini_length:

        print(d_slice.meta.endtime - d_slice.meta.starttime)
        print("SHORT")
        raise AttributeError

    p_slice = p.slice(d_slice.meta.starttime, d_slice.meta.endtime)


    # Must be a single trace
    if len(p_slice.traces) != 1:

        print("MUCH")
        raise AttributeError

    p_slice = p_slice.traces[0]

    sta = max(p_slice.meta.starttime, d_slice.meta.starttime)
    end = min(p_slice.meta.endtime, d_slice.meta.endtime)

# ################ modify if end - sta < 18000 - if end - sta < 6000 ################

    if end - sta < mini_length:

        print("SHORT")
        raise AttributeError

    return Stream().append(d_slice).slice(sta, end), Stream().append(p_slice).slice(sta, end)


def update_stats(tr, stla, stlo, stel, ws, ss):

    tr.stats.sac = AttribDict()
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    tr.stats.sac.user8 = ws
    tr.stats.sac.user9 = ss

    return tr


def load_db(fname):

    db = pickle.load(open(fname, 'rb'))
    # for k, v in db.items():

    #     db[k] = meta_data(v)

    return db


class meta_data(dict):

    def __init__(self, stats):

        self.__dict__ = self
        self.network = stats[0]
        self.station = stats[1]
        self.stla = stats[2]
        self.stlo = stats[3]
        self.stel = stats[4]
        self.azim = stats[5]
        self.cha = stats[6]
        self.dstart = stats[7]
        self.dend = stats[8]


if __name__ == '__main__':

    network = sys.argv[1]
    station = sys.argv[2]

    dbfile = '{}.{}.pkl'.format(network, station)
    stationdb = load_db(dbfile)

    process(stationdb, network, station)
