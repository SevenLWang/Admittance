import os
import fnmatch
from obspy import UTCDateTime
from obspy.core import read, Stream
import pickle
import numpy as np
import sys


def process(db, network, station):

    key = '{}.{}'.format(network, station)
    sta = db[key]  # sta is the meta_data

    UTC_sta = UTCDateTime("2018-01-01T00:00:00")
    UTC_end = UTCDateTime("2019-01-01T00:00:00")
    cur_time = UTC_sta

    ws = 2000      # window length
    dn = 1         # timestep
    frag_unit = 5  # number of segments

    freq = np.fft.fftfreq(ws, dn)

    days = open(sta.station + "/days", "w+")
    list_ratios = []

    while cur_time < UTC_end:

        julday = cur_time.julday
        year = cur_time.year
        trZ, trP = get_data(sta.station + '/Selected_Data/', year, julday)

        if len(trZ) < 1 or len(trP) < 1:

            cur_time += 3600 * 24
            continue

        else:

            ss = ws
            freq = np.fft.fftfreq(ws, dn)

            tx, nd_z = sliding_window(trZ[0].data, ws, ss)
            ty, nd_p = sliding_window(trP[0].data, ws, ss)
            nd = min(nd_z, nd_p)

            coh, ratio = window_selection_DPRatio(tx,ty,freq,nd,frag_unit)
            Ratio_Average = ratio_average(freq,ratio,ws,dn)

            if Ratio_Average:

                list_ratios.append(Ratio_Average)
                days.write(str(year)+"."+str(julday).zfill(3) + '\n')
                print("daily working on  " + str(year) + ' ' + str(julday))
                cur_time += 3600 * 24

            else:

                cur_time += 3600 * 24

    with open(sta.station+"/ratio", "w") as f:

        line = ""
        for i in range(len(list_ratios)):

            for j in range(len(list_ratios[i])):

                line += "{} ".format(list_ratios[i][j])
            
            line += "\n"

        f.write(line)


def ratio_average(freq, ratio, ws, dn):

    Ratio_Sum = []
    num = 0
    for i in range(len(ratio)):

        All_ratio = ratio[i]
        if str(All_ratio[0]) != 'nan':

            num += 1
            if i == 0:

                Ratio_Sum = All_ratio

            else:

                for j in range(len(Ratio_Sum)):

                    Ratio_Sum[j] += All_ratio[j]

    Ratio_Mean = []
    for k in range(len(Ratio_Sum)):

        Ratio_Mean.append(Ratio_Sum[k] / num)

    return Ratio_Mean


def coherence(Gxy, Gxx, Gyy):
    return np.abs(Gxy) ** 2 / (Gxx * Gyy)


def admittance(Gxy, Gxx):
    return np.abs(Gxy) / Gxx


def sliding_window(a, ws, ss=None):
    '''
    Parameters
        a  - a 1D array
        ws - the window size, in samples
        ss - the step size, in samples. If not provided, window and step size
             are equal.

    """ PA: This is not my function """

    '''

    if None is ss:
        # no step size was provided. Return non-overlapping windows
        ss = ws

    # Calculate the number of windows to return
    valid = len(a)  # how long the window can move Oct 29
    nd = (valid) // ws  # step size means how long one step can move Oct 29
    out = np.ndarray((nd, ws), dtype=a.dtype)


    for i in range(nd):  

        start = i * ws
        stop = start + ss
        out[i] = a[start: stop] * np.hanning(ws)

        # TODO change to Gaussian?
        # Return the Hanning window with the maximum value normalized to one Oct 29

    return out, nd  # 开始的时候会有一个taper，舍弃第一个窗户


def window_selection_DPRatio(segment_Z,segment_P,freq,nd,frag_unit):

    # 5 2000s windows as a fragment
    fragment_numbers = nd // frag_unit
    fragment_label = np.zeros((fragment_numbers,), dtype=bool)
    list_Ratio = []
    list_Coherence = []
    goodwins = 0
    for i in range(fragment_numbers-1):

        fragment_Z = segment_Z[i * frag_unit:(i+1) * frag_unit]
        fragment_P = segment_P[i * frag_unit:(i+1) * frag_unit]

        ftx = np.fft.fft(fragment_Z)
        fty = np.fft.fft(fragment_P)
        
        Gxx = np.abs(np.mean(np.conj(ftx) * ftx, axis=0))
        Gyy = np.abs(np.mean(np.conj(fty) * fty, axis=0))
        Gxy = np.mean(np.conj(ftx) * fty, axis=0)

        Coherence = coherence(Gxy, Gxx, Gyy)
        list_Coherence.append(Coherence)

        if Coherence[np.where(freq == 0.16)] >= 0.8 and Coherence[np.where(freq == 0.18)] >= 0.8 and Coherence[np.where(freq == 0.2)] >= 0.8 \
            and Coherence[np.where(freq == 0.22)] >= 0.8 and Coherence[np.where(freq == 0.17)] >= 0.8 and Coherence[np.where(freq == 0.19)] >= 0.8 \
                and Coherence[np.where(freq == 0.21)] >= 0.8:
        
            #### coherent band: 0.16Hz - 0.22Hz ####

                fragment_label[i] = True
                goodwins += 1

                uz = np.fft.fft(fragment_Z)
                p = np.fft.fft(fragment_P)

                numerator = uz * np.conj(p)
                denominator = p * np.conj(p)

                Ratio = np.abs(np.mean(numerator, axis=0) / np.mean(denominator, axis=0))
                list_Ratio.append(Ratio)

        else:

            fragment_label[i] = False

    return list_Coherence, list_Ratio


def get_data(filedir, year, julday):
    """
    Function to read all available receiver functions that meet SNR threshold
    :param filedir: String()
    :param julday: Integer
    :param year: Integer
    :return: Stream() objects
    """
    # Define empty streams
    trNZ = Stream()
    trNP = Stream()
    trZ = Stream()
    trP = Stream()

    filename = str(year) + "." + str(julday).zfill(3)

    # Loop through directory and load files
    for file in os.listdir(filedir):

        if fnmatch.fnmatch(file, filename + '.BHZ'):

            tr = read(filedir + file)
            trZ.append(tr[0])

        elif fnmatch.fnmatch(file, filename + '.P'):

            tr = read(filedir + file)
            trP.append(tr[0])

    return trZ, trP


# Loads station db and builds attribute dict of station stats
def load_db(fname):

    db = pickle.load(open(fname, 'rb'))

    # for k, v in db.items():  # k is station's name v is attributes of the station

    #     db[k] = meta_data(v)

    return db


# Attribute dict class
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



if __name__ == "__main__":

    network = sys.argv[1]
    station = sys.argv[2]

    dbfile = '{}.{}.pkl'.format(network, station)
    stationdb = load_db(dbfile)

    print("working on " + network + '.' + station)
    process(stationdb, network, station)

