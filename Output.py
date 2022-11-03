# coding=utf-8
from tkinter import NE
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import dates as mdates
import datetime
import copy
import sys
import pandas as pd
from tqdm import tqdm
from obspy import UTCDateTime
import datetime 
import os



def error_bar(ratios, freq, f, ws, dt):
    
    ste = []

    index = []
    for i in range(len(f)):

        start = np.where(freq == round(f[i] - 0.003, 3))[0]
        end = np.where(freq == round(f[i] + 0.004, 3))[0]
        index.append(np.arange(start, end))
        

    for i in range(len(index)):  # frequency

        values = []

        for r in ratios:  # day

            value = []

            for j in r[index[i]]:

                value.append(j)
            
            values.append(np.mean(value))

        std = np.std(values)  # standard deviation of each frequency in this week
        ste.append(std / np.sqrt(len(values)))  # stand error

    return ste


def plot_change_weekly(freq, ratios, days, ws, dt, title=None, size=15, sta=None):
    
    if sta == 'AXBA1':
        
        f = [0.14,0.15,0.16,0.17,0.18,0.19] # for AXBA1
        
    elif sta == 'AXCC1' or sta =='AXEC2':

        f = [0.16,0.17,0.18,0.19,0.20,0.21,0.22]   # for AXCC1 and AXEC2
        
    F = np.linspace(0.01, 0.3, 30)  # Frequnecy band during Calculation

    erp = 114  # julday of 2015.04.24 eruption 
    
    indexes = []
    week_ratios = []
    error_bars = [[] for i in range(len(f))]  # standard error
    Errorbars = [[] for i in range(len(F))]

    for i in range((365 * 5 + 366 * 2 + 70) // 7):  # Week numbers

        cur_index = []  # Index of the specific week
        cur_ratio = []
        start = i * 7 + 1
        end = start + 7

        if start < erp < end:  # Separately calculate the pre- and post-eruption admittance

            end = erp

        for j in range(start, end):

            try:
                
                x = days.index(j)
                cur_index.append(j)
                cur_ratio.append(ratios[x])

            except ValueError:  # None of this day

                pass

        if len(cur_index) != 0:

            mean = np.mean(cur_ratio[:], axis=0)
            week_ratios.append(mean)
            indexes.append(start)

            ste = error_bar(cur_ratio[:], freq, f, ws, dt)  
            Ste = error_bar(cur_ratio[:], freq, F, ws, dt)

            for k in range(len(error_bars)):

                error_bars[k].append(ste[k])
            
            for k in range(len(Errorbars)):
                
                Errorbars[k].append(Ste[k])


        if end == erp:  # Special for the post-eruption week

            cur_index = []
            cur_ratio = []

            for j in range(end, start + 7):

                try:

                    x = days.index(j)
                    cur_index.append(j)
                    cur_ratio.append(ratios[x])

                except ValueError:

                    pass

            if len(cur_index) != 0:

                mean = np.mean(cur_ratio[:], axis=0)
                week_ratios.append(mean)
                indexes.append(end)
                ste = error_bar(cur_ratio[:], freq, f, ws, dt)
                Ste = error_bar(cur_ratio[:], freq, F, ws, dt)
                
                for k in range(len(error_bars)):

                    error_bars[k].append(ste[k])
                
                for k in range(len(Errorbars)):

                    Errorbars[k].append(Ste[k])


    return plot_time_change(freq, week_ratios[:], indexes, ws, dt, sta, title, size, error_bars, Errorbars, f, F)


def plot_time_change(freq, ratios, days, ws, dt, key=None, title=None, size=15, err_bars=None, Err_bars=None, f=[], F=[]):

    date_start = datetime.datetime(2015,1,1,0,0,0)
    delta = datetime.timedelta(days=1)
    dates = []

    for d in days:

        d = int(d)
        cur_day = date_start + delta*(d-1)
        dates.append(cur_day)

    if f == []:

        if key == 'AXBA1':

            f = [0.14,0.15,0.16,0.17,0.18,0.19]

        elif key == 'AXCC1' or key =='AXEC2':

            f = [0.16,0.17,0.18,0.19,0.20,0.21,0.22]


    index = []
    for i in range(len(f)):

        start = np.where(freq == round(f[i] - 0.003, 3))[0]
        end = np.where(freq == round(f[i] + 0.004, 3))[0]
        index.append(np.arange(start, end))



    ### Output ###
    
    os.system("rm -r {}/{}_*Hz_ratio".format(key, key))
    F = np.linspace(0.01, 0.3, 30)
    AVG = [[] for i in range(len(F))]
    Errbars = [[] for i in range(len(F))]

    Index = []
    for i in range(len(F)):

        start = np.where(freq == round(F[i] - 0.003, 3))[0]
        end = np.where(freq == round(F[i] + 0.004, 3))[0]
        Index.append(np.arange(start, end))


    week_days = open("{}/{}_week_days".format(key, key),"w")
    for j in range(len(ratios)):

        if Err_bars[0][j] != 0:

            d = int(days[j])
            cur_day = date_start + delta*(d-1)
            week_days.write("{}\n".format(str(cur_day)))

            for i in range(len(F)):
                    
                value = np.average(ratios[j][Index[i]])
                AVG[i].append(value)

                value = Err_bars[i][j]
                Errbars[i].append(value)


    np.savetxt("{}/{}_week_ratio".format(key, key),AVG)
    np.savetxt("{}/{}_week_errorbar".format(key, key),Errbars)

    for i in range(len(F)):

        f_ratio = open("%s/%s_%.2fHz_ratio" % (key, key, F[i]), "w")

        for j in range(len(AVG[0])):

            date = dates[j]
            date = date.strftime("%Y-%m-%dT%H:%M:%S")

            line = "{} {} {}\n".format(date, AVG[i][j], Err_bars[i][j])
            f_ratio.write(line)
    
    
    
    ### Plot ###

    fig = plt.figure(dpi=110)
    labels = []
    colors = ['b', 'g', 'y', 'c', 'm', 'r', 'tan', 'k']

    avg = [[] for i in range(len(f))]
    for r in ratios:

        for i in range(len(f)):

            value = np.average(r[index[i]])
            avg[i].append(value)    # week average admittance

    for i in range(len(f)):

        if err_bars:

            plt.errorbar(dates, avg[i], yerr=err_bars[i], fmt='o', markersize=2, elinewidth=1, capsize=2, capthick=1)
            
        else:

            plt.plot(dates, np.log10(avg[i]), linestyle='-', marker='o',markersize=3)
            
        labels.append(f[i])

    plt.ylim()
    plt.grid(ls='--')

    ax = plt.gca()
    locator = mdates.AutoDateLocator()
    formatter = mdates.DateFormatter('%Y-%b-%d')
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    fig.autofmt_xdate()
    
    font = {'weight' : 'normal', 'fontsize' : 'large'}
    plt.xlabel("Date",font)
    plt.ylabel("D/P-Ratio,  $\mathregular{log_{10}}$ (m / GPa)",font)
    

    plt.tight_layout()
    name = '%s.png' % (title)
    plt.savefig(name,bbox_inches = 'tight')

    plt.show()

    return fig


if __name__ == "__main__":

    network = sys.argv[1]
    station = sys.argv[2]

    file_days = open("%s/days" % (station))
    days = []
    for i in file_days:
        year = int(i[0:4])
        day = int(i[5:-1])

        if 2016 < year < 2021:
            days.append(731 + 365*(year-2017) + day)
        elif year <= 2016:
            days.append(365*(year-2015) + day)
        elif year >= 2021:
            days.append(2192 + 365*(year-2021) + day)

    ratios = np.loadtxt("%s/ratio" % (station))
    ws = len(ratios[0])
    dn = 1   # time step
    multiple = ws * dn
    freq = np.fft.fftfreq(ws, dn)

    for i in range(len(freq)):

        freq[i] = round(freq[i], 5)

    title = '%s_%s_Weekly_2015-2022' % (network, station)

    plot_change_weekly(freq, ratios, days, ws, dn, title=title, sta=station)

