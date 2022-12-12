# Displacement - Pressure Ratio

## Basic Pre-Process

- Remove mean
- Remove trend
- Taper
- Remove instrument response
- Filt & Decimation

## DP-Ratio Calculation

- Query for the station pickle file (need obstools python package)


- Data filteration (Put the original data in the $station/Raw_Data, and create the $station/Selected_Data)

    - Set the starttime, endtime and minimal data length in Cut_DailyData.py
    - python Cut_DailyData.py $network $station


- Ratio calculation

    - Set the starttime and endtime in Calculate_Ratio.py
    - Set the window length, sampling rate and fragment number
    - python Calculate_Ratio.py $network $station

- Ratio output

    - python Output.py $network $station





## 