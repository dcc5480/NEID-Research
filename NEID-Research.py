import csv
import matplotlib.pyplot as plt
import math
import random
import numpy as np
import time
import os
#https://github.com/tamarervin/SolAster #https://tamarervin.github.io/SolAster/about/
import astropy
import SolAster
import SolAster.tools.calculation_funcs as sfuncs
import SolAster.tools.lbc_funcs as lbfuncs
import SolAster.tools.coord_funcs as ctfuncs
import sunpy
from sunpy.net import Fido, attrs as a
import drms
from datetime import datetime,timedelta,timezone
import datetime
import pandas as pd
import sunpy.map
from sunpy.coordinates import frames
import sys
sys.path.append(os.path.realpath('../../'))
import SolAster.tools.rvs as rvs
import SolAster.tools.utilities as utils
from SolAster.tools.settings import *
from SolAster.tools.plotting_funcs import hmi_plot
import scipy
from scipy.signal import savgol_filter


def install(module):import subprocess;import sys;exec(f'subprocess.check_call([sys.executable, "-m", "pip", "install", "{module}"])')

with open("/Users/suoenallecsim/Documents/combined_rvs_espressoG2_norm=cont&mask=1.5_fwtf=1.0.csv") as file:
  L=csv.reader(file)
  L=[next(L)]+[[i[j]if not i[j] or j in(6,18,20,21,26,43,46,62,63,79) else float(i[j])for j in range(len(i))] for i in L]
  Lx=[L[i][22]for i in range(1,len(L))]
  Ly=[L[i][23]for i in range(1,len(L))]
  Lca1=[L[i][4]for i in range(1,len(L))]
  Lca2=[L[i][5]for i in range(1,len(L))]
  Lha01=[L[i][7]for i in range(1,len(L))]
  Lha02=[L[i][8]for i in range(1,len(L))]
  Lha11=[L[i][9]for i in range(1,len(L))]
  Lha12=[L[i][10]for i in range(1,len(L))]
  Lhe1=[L[i][11]for i in range(1,len(L))]
  Lhe2=[L[i][12]for i in range(1,len(L))]
  Lna=[L[i][13]for i in range(1,len(L))]

for i in range(1):
  Lx1=[]
  Ly1=[]
  Lca11=[]
  Lha011=[]
  Lha121=[]
  Lhe11=[]
  t1=int(Lx[0])+1
  t2=0
  for i in range(len(Lx)):
    if Lx[i]>=t1:
      Lx1.append(Lx[t2:i])
      Ly1.append(Ly[t2:i])
      Lca11.append(Lca1[t2:i])
      Lha011.append(Lha01[t2:i])
      Lha121.append(Lha12[t2:i])
      Lhe11.append(Lhe1[t2:i])
      t1=int(Lx[i])+1
      t2=i
  Lx1.append(Lx[t2:])
  Ly1.append(Ly[t2:])
  Lca11.append(Lca1[t2:])
  Lha011.append(Lha01[t2:])
  Lha121.append(Lha12[t2:])
  Lhe11.append(Lhe1[t2:])
  
  t1=int(len(Lx)/len(Lx1)/3)
  
  Lx2=[i for i in Lx1 if len(i)>t1]
  Ly2=[i for i in Ly1 if len(i)>t1]
  Lca12=[i for i in Lca11 if len(i)>t1]
  Lha012=[i for i in Lha011 if len(i)>t1]
  Lha122=[i for i in Lha121 if len(i)>t1]
  Lhe12=[i for i in Lhe11 if len(i)>t1]

  Lx2i=[int(i[0])for i in Lx2]
  Ly2a=[sum(i)/len(i)for i in Ly2]
  Lx3=[j for i in Lx2 for j in i]
  Ly3=[j for i in Ly2 for j in i]
  Lca13=[j for i in Lca12 for j in i]
  Lha013=[j for i in Lha012 for j in i]
  Lha123=[j for i in Lha122 for j in i]
  Lhe13=[j for i in Lhe12 for j in i]
  

#plot side by side with separate subplots
def multiPlot(L3,L2='default',L1='default',vert=True):
  if L1=='default':L1=Lx
  if L2=='default':L2=Ly
  if vert:t1,t2=2,1
  else:t1,t2=1,2
  plt.subplot(t1,t2,1)
  plt.plot(L1,L2)
  plt.ylabel("L2")
  plt.subplot(t1,t2,2)
  plt.plot(L1,L3)
  plt.ylabel("L3")
  plt.show()
  #generalize this
  #Lma=movingAverage(20,Ly);Lha01ma=movingAverage(20,Lha01);plt.subplot(2,1,1);plt.plot(Lx,Ly);plt.plot(Lx,Lma);plt.subplot(2,1,2);plt.plot(Lx,Lha01);plt.plot(Lx,Lha01ma);plt.show()

def movingAverage(n,L2):
  #uses consecutive groups of n data points, disregarding time
  #temp hard-coded value:
  n=10
  
  L3=L2[:n//2]
  for i in range(n,len(L2)):
    L3.append(sum(L2[j] for j in range(i-n,i))/n)
  L3+=L2[-n+n//2:]
  return L3

def derivative(L):
  L1=[L[0]]
  for i in range(1,len(L)):
    L1.append((L[i]-L[i-1])/(Lx[i]-Lx[i-1]))
  return L1

def zeroCount(L,days=Lx1):
  #counts zeros, can be used with derivative to count peaks
  L1=[]
  c=0
  s=1
  for i in days:
    t1=0
    for j in range(len(i)):
      if s*L[c+j]<=0:
        t1+=1
        s=-s
    L1.append(t1)
    c+=len(i)
  return L1

def closestListValue(x,L,low=0,hig=None):
  #requires sorted L
  #binary search
  if not low:low=0
  else:low=max(0,low)
  if not hig:hig=len(L)-1
  else:hig=min(hig,len(L)-1)
  while abs(low-hig)>1:
    mid=(low+hig)//2
    if x<L[mid]:hig=mid
    else:low=mid
  if abs(x-L[low])<=abs(x-L[hig]):return low
  return hig

def present1():
  plt.rcParams.update({'font.size': 14})
  #raw data with 20 point double moving averages (ma) overlayed on left
  #20 point double ma of derivative of 20 point ma on right with 0 overlayed
  #top is rv vs jd and bottom is ca1 vs jd
  ax1=plt.subplot(2,2,1);plt.plot(Lx,Ly);plt.plot(Lx,Lymama);plt.xlabel("JD (days)");plt.ylabel("RV (m/s) with 10pt double moving average overlayed");plt.subplot(2,2,2,sharex=ax1);plt.plot(Lx,Lymadmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.xlabel("JD (days)");plt.ylabel("∆RV (m/s/day) smoothed and with 0 overlayed");plt.subplot(2,2,3,sharex=ax1);plt.plot(Lx,Lca1);plt.plot(Lx,Lca1mama);plt.xlabel("JD (days)");plt.ylabel("RV (m/s) of Ca1, smoothed");plt.subplot(2,2,4,sharex=ax1);plt.plot(Lx,Lca1madmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.xlabel("JD (days)");plt.ylabel("∆RV (m/s/day) of Ca1, smoothed");plt.show()
  #same but with ha01 instead of ca1
  ax1=plt.subplot(2,2,1);plt.plot(Lx,Ly);plt.plot(Lx,Lymama);plt.xlabel("JD (days)");plt.ylabel("RV (m/s) with 10pt double moving average overlayed");plt.subplot(2,2,2,sharex=ax1);plt.plot(Lx,Lymadmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.xlabel("JD (days)");plt.ylabel("∆RV (m/s/day) smoothed and with 0 overlayed");plt.subplot(2,2,3,sharex=ax1);plt.plot(Lx,Lha01);plt.plot(Lx,Lha01mama);plt.xlabel("JD (days)");plt.ylabel("RV (m/s) of Ha01, smoothed");plt.subplot(2,2,4,sharex=ax1);plt.plot(Lx,Lha01madmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.xlabel("JD (days)");plt.ylabel("∆RV (m/s/day) of Ha01, smoothed");plt.show()
  #same but with counts of zeros rather than the derivative
  #ax1=plt.subplot(2,2,1);plt.plot(Lx,Ly);plt.plot(Lx,Lymama);plt.subplot(2,2,2,sharex=ax1);plt.plot(Lx11,Lymadmamaz);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.subplot(2,2,3,sharex=ax1);plt.plot(Lx,Lca1);plt.plot(Lx,Lca1mama);plt.subplot(2,2,4,sharex=ax1);plt.plot(Lx11,Lca1madmamaz);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.show()
  #histogram of counts of zeros per day for rv and for ca1
  ax1=plt.subplot(1,2,1);plt.hist(Lymadmamazs[0],bins=len(Lymadmamazs[0]),weights=Lymadmamazs[1]);plt.xlabel("Count of Zeros in smoothed ∆RV per Day");plt.ylabel("Count of Days");plt.subplot(1,2,2,sharex=ax1);plt.hist(Lca1madmamazs[0],bins=len(Lca1madmamazs[0]),weights=Lca1madmamazs[1]);plt.xlabel("Count of Zeros in smoothed ∆RV for Ca1 per Day");plt.ylabel("Count of Days");plt.show()
  #same but with ha01
  ax1=plt.subplot(1,2,1);plt.hist(Lymadmamazs[0],bins=len(Lymadmamazs[0]),weights=Lymadmamazs[1]);plt.xlabel("Count of Zeros in smoothed ∆RV per Day");plt.ylabel("Count of Days");plt.subplot(1,2,2,sharex=ax1);plt.hist(Lha01madmamazs[0],bins=len(Lha01madmamazs[0]),weights=Lha01madmamazs[1]);plt.xlabel("Count of Zeros in smoothed ∆RV for Ha01 per Day");plt.ylabel("Count of Days");plt.show()
  #heatmap of counts of zeros per day for ca1 vs rv and ha01 vs rv
  ax1=plt.subplot(1,2,1);plt.hist2d(Lymadmamaz,Lca1madmamaz,bins=(max(Lymadmamaz)-min(Lymadmamaz),max(Lca1madmamaz)-min(Lca1madmamaz)));plt.xlabel("Zeros in smoothed ∆RV");plt.ylabel("Zeros in smoothed ∆RV for Ca1");plt.subplot(1,2,2,sharey=ax1);plt.hist2d(Lymadmamaz,Lha01madmamaz,bins=(max(Lymadmamaz)-min(Lymadmamaz),max(Lha01madmamaz)-min(Lha01madmamaz)));plt.xlabel("Zeros in smoothed ∆RV");plt.ylabel("Zeros in smoothed ∆RV for Ha01");plt.show()
  #same but with low data point days removed (>(average per day)/3=52)
  #and bottom same but ha12 vs rv and he1 vs rv
  ax1=plt.subplot(2,2,1);plt.hist2d(Ly3madmamaz,Lca13madmamaz,bins=(max(Ly3madmamaz)-min(Ly3madmamaz),max(Lca13madmamaz)-min(Lca13madmamaz)));plt.xlabel("Zeros in smoothed ∆RV (no bad)");plt.ylabel("Zeros in smoothed ∆RV for Ca1 (no bad)");plt.subplot(2,2,2,sharey=ax1);plt.hist2d(Ly3madmamaz,Lha013madmamaz,bins=(max(Ly3madmamaz)-min(Ly3madmamaz),max(Lha013madmamaz)-min(Lha013madmamaz)));plt.xlabel("Zeros in smoothed ∆RV (no bad)");plt.ylabel("Zeros in smoothed ∆RV for Ha01 (no bad)");plt.subplot(2,2,3,sharey=ax1);plt.hist2d(Ly3madmamaz,Lha123madmamaz,bins=(max(Ly3madmamaz)-min(Ly3madmamaz),max(Lha123madmamaz)-min(Lha123madmamaz)));plt.xlabel("Zeros in smoothed ∆RV (no bad)");plt.ylabel("Zeros in smoothed ∆RV for Ha12 (no bad)");plt.subplot(2,2,4,sharey=ax1);plt.hist2d(Ly3madmamaz,Lhe13madmamaz,bins=(max(Ly3madmamaz)-min(Ly3madmamaz),max(Lhe13madmamaz)-min(Lhe13madmamaz)));plt.xlabel("Zeros in smoothed ∆RV (no bad)");plt.ylabel("Zeros in smoothed ∆RV for He1 (no bad)");plt.show()
  plt.rcParams.update({'font.size': 10})
  
def dataDownload(jd):
  #low,hig,mid:
  #Ly3madmamaz[16],Ly3madmamaz[442],Ly3madmamaz[227]
  t1=(datetime(1858,11,17,tzinfo=timezone.utc)+timedelta(jd-2400000.5)).strftime("%Y-%m-%d")
  Fido.fetch(Fido.search(a.Time(t1+'T00:00:00',t1+'T00:00:00')|a.Time(t1+'T12:00:00',t1+'T12:00:00')|a.Time(t1+'T23:59:59',t1+'T23:59:59'),a.jsoc.Series('hmi.ic_45s'),a.jsoc.Segment('continuum'),a.jsoc.Notify('dcc5480@psu.edu')))
  Fido.fetch(Fido.search(a.Time(t1+'T00:00:00',t1+'T00:00:00')|a.Time(t1+'T12:00:00',t1+'T12:00:00')|a.Time(t1+'T23:59:59',t1+'T23:59:59'),a.jsoc.Series('hmi.v_45s'),a.jsoc.Segment('Dopplergram'),a.jsoc.Notify('dcc5480@psu.edu')))
  Fido.fetch(Fido.search(a.Time(t1+'T00:00:00',t1+'T00:00:00')|a.Time(t1+'T12:00:00',t1+'T12:00:00')|a.Time(t1+'T23:59:59',t1+'T23:59:59'),a.jsoc.Series('hmi.m_45s'),a.jsoc.Segment('magnetogram'),a.jsoc.Notify('dcc5480@psu.edu')))


  
  

#noteworty areas:
  #2459930 Lca1,ca2,he1,...
  #2460030 to the end
  #2459710 to 2459750
  #2459380 to 2459500
  #2460043 to 2460045 using ma and ca1ma
  #2459567 to 2459569 using ma and ca1ma
  #2459269 to 2459273 using ma and ca1ma
  #2459622 to 2459624 using ma and ha01ma, note that rv and ha01 seem inversely correlated
#moving average of 20 seems effective
  #comparing ma of rv with ca1 shows possible peaks in ca1 <30 min before peaks in rv

#Ca S index

for i in range(1):
  #Lymad=derivative(Lyma);"etc"
  Lyma=movingAverage(20,Ly);Lymama=movingAverage(20,Lyma);Lymad=derivative(Lyma);Lymadma=movingAverage(20,Lymad);Lymadmama=movingAverage(20,Lymadma)
  Lca1ma=movingAverage(20,Lca1);Lca1mama=movingAverage(20,Lca1ma);Lca1mad=derivative(Lca1ma);Lca1madma=movingAverage(20,Lca1mad);Lca1madmama=movingAverage(20,Lca1madma)
  Lha01ma=movingAverage(20,Lha01);Lha01mama=movingAverage(20,Lha01ma);Lha01mad=derivative(Lha01ma);Lha01madma=movingAverage(20,Lha01mad);Lha01madmama=movingAverage(20,Lha01madma)
  #ax1=plt.subplot(2,2,1);plt.plot(Lx,Ly);plt.plot(Lx,Lymama);plt.subplot(2,2,2,sharex=ax1);plt.plot(Lx,Lymadmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.subplot(2,2,3,sharex=ax1);plt.plot(Lx,Lca1);plt.plot(Lx,Lca1mama);plt.subplot(2,2,4,sharex=ax1);plt.plot(Lx,Lca1madmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.show()
  #ax1=plt.subplot(2,2,1);plt.plot(Lx,Ly);plt.plot(Lx,Lymama);plt.subplot(2,2,2,sharex=ax1);plt.plot(Lx,Lymadmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.subplot(2,2,3,sharex=ax1);plt.plot(Lx,Lha01);plt.plot(Lx,Lha01mama);plt.subplot(2,2,4,sharex=ax1);plt.plot(Lx,Lha01madmama);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.show()

  #ax1=plt.subplot(2,2,1);plt.plot(Lx,Ly);plt.plot(Lx,Lymama);plt.subplot(2,2,2,sharex=ax1);plt.plot(Lx11,Lymadmamaz);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.subplot(2,2,3,sharex=ax1);plt.plot(Lx,Lca1);plt.plot(Lx,Lca1mama);plt.subplot(2,2,4,sharex=ax1);plt.plot(Lx11,Lca1madmamaz);plt.plot([Lx[0],Lx[-1]],[0,0]);plt.show()

  Lx11=[int(i[0])for i in Lx1]
  Lymadmamaz=zeroCount(Lymadmama)
  Lca1madmamaz=zeroCount(Lca1madmama)
  Lha01madmamaz=zeroCount(Lha01madmama)
  Lymadmamazs=list(zip(*[(i,Lymadmamaz.count(i))for i in set(Lymadmamaz)]))
  Lca1madmamazs=list(zip(*[(i,Lca1madmamaz.count(i))for i in set(Lca1madmamaz)]))
  Lha01madmamazs=list(zip(*[(i,Lha01madmamaz.count(i))for i in set(Lha01madmamaz)]))
  #ax1=plt.subplot(1,2,1);plt.hist(Lymadmamazs[0],bins=len(Lymadmamazs[0]),weights=Lymadmamazs[1]);plt.subplot(1,2,2,sharex=ax1);plt.hist(Lca1madmamazs[0],bins=len(Lca1madmamazs[0]),weights=Lca1madmamazs[1]);plt.show()
  #ax1=plt.subplot(1,2,1);plt.hist(Lymadmamazs[0],bins=len(Lymadmamazs[0]),weights=Lymadmamazs[1]);plt.subplot(1,2,2,sharex=ax1);plt.hist(Lha01madmamazs[0],bins=len(Lha01madmamazs[0]),weights=Lha01madmamazs[1]);plt.show()

  Ly3ma=movingAverage(20,Ly3);Ly3mama=movingAverage(20,Ly3ma);Ly3mad=derivative(Ly3ma);Ly3madma=movingAverage(20,Ly3mad);Ly3madmama=movingAverage(20,Ly3madma)
  Lca13ma=movingAverage(20,Lca13);Lca13mama=movingAverage(20,Lca13ma);Lca13mad=derivative(Lca13ma);Lca13madma=movingAverage(20,Lca13mad);Lca13madmama=movingAverage(20,Lca13madma)
  Lha013ma=movingAverage(20,Lha013);Lha013mama=movingAverage(20,Lha013ma);Lha013mad=derivative(Lha013ma);Lha013madma=movingAverage(20,Lha013mad);Lha013madmama=movingAverage(20,Lha013madma)
  Lha123ma=movingAverage(20,Lha123);Lha123mama=movingAverage(20,Lha123ma);Lha123mad=derivative(Lha123ma);Lha123madma=movingAverage(20,Lha123mad);Lha123madmama=movingAverage(20,Lha123madma)
  Lhe13ma=movingAverage(20,Lhe13);Lhe13mama=movingAverage(20,Lhe13ma);Lhe13mad=derivative(Lhe13ma);Lhe13madma=movingAverage(20,Lhe13mad);Lhe13madmama=movingAverage(20,Lhe13madma)


  Ly3madmamaz=zeroCount(Ly3madmama,Lx2)
  Lca13madmamaz=zeroCount(Lca13madmama,Lx2)
  Lha013madmamaz=zeroCount(Lha013madmama,Lx2)
  Lha123madmamaz=zeroCount(Lha123madmama,Lx2)
  Lhe13madmamaz=zeroCount(Lhe13madmama,Lx2)



















def fullPipeline(start_date,end_date,jd=False,cadence=24*60*60//3,name="example"):
  #input start_date and end_date as eg. datetime.datetime(2021, 1, 1, 0, 0, 0)
  #https://github.com/tamarervin/SolAster/blob/main/SolAster/examples/full_pipeline.ipynb
  class Inputs:
    """
    Class to hold user specified inputs to run examples.
    See README or documentation site for additional information.
    """
    randname=random.randint(0,10**12)
    # name of csv file to store calculations
    csv_name = f'{name}_calcs_{randname}.csv'

    # name of instrument to use for calculation of RV model
    # choose either 'NEID' or 'HARPS-N'
    inst = 'NEID'

    # querying cadence in seconds
    #cadence = 24 * 60 * 60//3

    # start date for calculationsx
    #start_date = datetime.datetime(2021, 1, 1, 0, 0, 0)

    # end date for calculations
    #end_date = datetime.datetime(2021, 1, 2, 0, 0, 0)

    # True if outputting diagnostic plots
    diagnostic_plots = True
    # path to save diagnostic figure or none
    save_fig = None

    # figure title for plot of calculations
    fig_title = f'{name}_plot_{randname}.png'
  if jd:
    start_date=datetime.datetime(1858,11,17)+timedelta(start_date-2400000.5)
    end_date=datetime.datetime(1858,11,17)+timedelta(end_date-2400000.5)

  for i in '1':
    # check input formats

    start_date, end_date, cadence, csv_name = utils.check_inputs(CsvDir.CALC,  start_date,  end_date, cadence,  Inputs.csv_name)


    # print out csv title
    print("Beginning calculation of values for csv file: " + csv_name)

    # List of header strings
    row_contents = ['date_obs', 'date_jd', 'rv_model', 'v_quiet', 'v_disc', 'v_phot', 'v_conv', 'f_bright', 'f_spot', 'f',
        'Bobs', 'vphot_bright', 'vphot_spot', 'f_small', 'f_large', 'f_network', 'f_plage',
        'quiet_flux', 'ar_flux', 'conv_flux', 'pol_flux', 'pol_conv_flux', 'vconv_quiet', 'vconv_large',
        'vconv_small']

    # create file names
    csv_file = os.path.join(CsvDir.CALC, csv_name)
    bad_dates_csv = os.path.join(CsvDir.CALC, csv_name[:-4]+'_bad_dates.csv')
    print(bad_dates_csv)
    utils.append_list_as_row(csv_file, row_contents)


    # get hmi data products
    time_range = datetime.timedelta(seconds=22)
    physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

    # get dates list
    xy = (end_date - start_date).seconds + (end_date - start_date).days * 24 * 3600
    dates_list = [start_date + datetime.timedelta(seconds=cadence*x) for x in range(0, int(xy/cadence))]
  exec("""for i, date in enumerate(dates_list):
    # convert the date to a string -- required for use in csv file
    date_str, date_obj, date_jd = utils.get_dates(date)

    # pull image within specified time range
    result = Fido.search(a.Time(str(date_obj - time_range), str(date_obj + time_range)),
                         a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2])

    # add file to list
    file_download = Fido.fetch(result)

    # remove unusable file types
    good_files = []
    for file in file_download:
        name, extension = os.path.splitext(file)
        if extension == '.fits':
            good_files.append(file)

    if len(good_files) != 3:
        # add the data
        # append these values to the csv file
        save_vals = [date_str, 'not three good files']
        # print that the files are missing
        print('Not three good files: ' + date_str + ' index: ' + str(i))

        pass
    else:
        # convert to map sequence
        map_seq = sunpy.map.Map(sorted(good_files))

        # check for missing data types
        missing_map = False
        # split into data types
        for j, map_obj in enumerate(map_seq):
            if map_obj.meta['content'] == 'DOPPLERGRAM':
                vmap = map_obj
            elif map_obj.meta['content'] == 'MAGNETOGRAM':
                mmap = map_obj
            elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
                imap = map_obj
            else:
                missing_map = True

        if missing_map:
            print("Missing a data product for " + date_str)

            # add the data
            # append these values to the csv file
            save_vals = [date_str, 'missing data product']
            pass

        else:
            # coordinate transformation for maps
            x, y, pdim, r, d, mu = ctfuncs.coordinates(vmap)
            wij, nij, rij = ctfuncs.vel_coords(x, y, pdim, r, vmap)

            # remove bad mu values
            vmap, mmap, imap = ctfuncs.fix_mu(mu, [vmap, mmap, imap], mu_cutoff=Parameters.mu_cutoff)

            # calculate relative positions
            deltaw, deltan, deltar, dij = sfuncs.rel_positions(wij, nij, rij, vmap)

            # calculate spacecraft velocity
            vsc = sfuncs.spacecraft_vel(deltaw, deltan, deltar, dij, vmap)

            # optimized solar rotation parameters
            a_parameters = [Parameters.a1, Parameters.a2, Parameters.a3]

            # calculation of solar rotation velocity
            vrot = sfuncs.solar_rot_vel(wij, nij, rij, deltaw, deltan, deltar, dij, vmap, a_parameters)

            # calculate corrected velocity
            corrected_vel = vmap.data - np.real(vsc) - np.real(vrot)

            # corrected velocity maps
            map_vel_cor = sfuncs.corrected_map(corrected_vel, vmap, map_type='Corrected-Dopplergram',
                                               frame=frames.HeliographicCarrington)

            # limb brightening
            Lij = lbfuncs.limb_polynomial(imap)

            # calculate corrected data
            Iflat = imap.data / Lij

            # corrected intensity maps
            map_int_cor = sfuncs.corrected_map(Iflat, imap, map_type='Corrected-Intensitygram',
                                               frame=frames.HeliographicCarrington)

            # calculate unsigned field strength
            Bobs, Br = sfuncs.mag_field(mu, mmap, B_noise=Parameters.B_noise, mu_cutoff=Parameters.mu_cutoff)

            # corrected observed magnetic data map
            map_mag_obs = sfuncs.corrected_map(Bobs, mmap, map_type='Corrected-Magnetogram',
                                               frame=frames.HeliographicCarrington)

            # radial magnetic data map
            map_mag_cor = sfuncs.corrected_map(Br, mmap, map_type='Corrected-Magnetogram',
                                               frame=frames.HeliographicCarrington)

            # calculate magnetic threshold
            active, quiet = sfuncs.mag_thresh(mu, mmap, Br_cutoff=Parameters.Br_cutoff, mu_cutoff=Parameters.mu_cutoff)

            # calculate intensity threshold
            fac_inds, spot_inds = sfuncs.int_thresh(map_int_cor, active, quiet)

            # create threshold array
            thresh_arr = sfuncs.thresh_map(fac_inds, spot_inds)

            # full threshold maps
            map_full_thresh = sfuncs.corrected_map(thresh_arr, mmap, map_type='Threshold',
                                                   frame=frames.HeliographicCarrington)

            # create diagnostic plots
            if i == 0 and Inputs.diagnostic_plots == True:
                    hmi_plot(map_int_cor, map_mag_cor, map_vel_cor, fac_inds, spot_inds, mu, save_fig=Inputs.save_fig)

            ### velocity contribution due to convective motion of quiet-Sun
            v_quiet = sfuncs.v_quiet(map_vel_cor, imap, quiet)

            ### velocity contribution due to rotational Doppler imbalance of active regions (faculae/sunspots)
            # calculate photospheric velocity
            v_phot, vphot_bright, vphot_spot = sfuncs.v_phot(quiet, active, Lij, vrot, imap, mu, fac_inds, spot_inds, mu_cutoff=Parameters.mu_cutoff)

            ### velocity contribution due to suppression of convective blueshift by active regions
            # calculate disc-averaged velocity
            v_disc = sfuncs.v_disc(map_vel_cor, imap)

            # calculate convective velocity
            v_conv = v_disc - v_quiet

            ### filling factor
            # calculate filling factor
            f_bright, f_spot, f = sfuncs.filling_factor(mu, mmap, active, fac_inds, spot_inds, mu_cutoff=Parameters.mu_cutoff)

            ### unsigned magnetic flux
            # unsigned observed flux
            unsigned_obs_flux = sfuncs.unsigned_flux(map_mag_obs, imap)

            ### calculate the area filling factor
            pixA_hem = ctfuncs.pix_area_hem(wij, nij, rij, vmap)
            area = sfuncs.area_calc(active, pixA_hem)
            f_small, f_large, f_network, f_plage = sfuncs.area_filling_factor(active, area, mu, mmap, fac_inds,
                                                                              athresh=Parameters.athresh,
                                                                              mu_cutoff=Parameters.mu_cutoff)

            ### get the unsigned flux
            quiet_flux, ar_flux, conv_flux, pol_flux, pol_conv_flux = sfuncs.area_unsigned_flux(map_mag_obs, imap,
                                                                                                    area,
                                                                                                    active,
                                                                                                athresh=Parameters.athresh)

            ### get area weighted convective velocities
            vconv_quiet, vconv_large, vconv_small = sfuncs.area_vconv(map_vel_cor, imap, active, area, athresh=Parameters.athresh)
            ### calculate model RV
            rv_model = rvs.calc_model(Inputs.inst, v_conv, v_phot)

            # intensity flux to check
            int_flux = np.nansum(imap.data)

            # make array of what we want to save
            save_vals = [rv_model, v_quiet, v_disc, v_phot, v_conv, f_bright, f_spot, f, unsigned_obs_flux, vphot_bright,
                         vphot_spot, f_small, f_large, f_network, f_plage, quiet_flux, ar_flux,
                         conv_flux, pol_flux, pol_conv_flux, vconv_quiet, vconv_large, vconv_small, int_flux]

            # round stuff
            rounded = np.around(save_vals, 3)
            round_vals = [date_str, date_jd]
            for val in rounded:
                round_vals.append(val)

            # append these values to the csv file
            utils.append_list_as_row(csv_file, round_vals)

            # print that the date is completed
            print('Calculations and save to file complete for ' + date_str + ' index: ' + str(i))""")
  csv_file = os.path.join(CsvDir.CALC, Inputs.csv_name)
  for i in '1':
      component_df = pd.read_csv(csv_file)

      date_jd = component_df.date_jd.values
      x = date_jd - date_jd[0]
      y_list = [component_df.f.values, component_df.Bobs.values, component_df.v_conv.values, component_df.v_phot.values,
          component_df.rv_model.values - np.median(component_df.rv_model.values)]

      # plot labels
      xlabel = 'Days since ' + str(int(date_jd[0])) + ' JD'
      ylabel_list = [r'$\rm f$' '\n' r'$\rm$[%]',
               r'$\rm B_{\rm obs}$' '\n' r'$\rm [G]$',
               r'$\rm v_{\rm conv}$' '\n' r'$\rm[m s^{-1}]$',
               r'$\rm v_{\rm phot}$' '\n' r'$\rm[m s^{-1}]$',
               r'$\rm RV_{\rm model}$' '\n' r'$\rm[m s^{-1}]$']

      # set up figure
      fig, axs = plt.subplots(len(y_list), 1, sharex='all', figsize=[6, 1.5 * len(y_list)],  gridspec_kw={'hspace': 0})

      # set up axes labels
      for i in range(0, len(axs)):
        axs[i].set(ylabel=ylabel_list[i])
        rng = (y_list[i].max() - y_list[i].min())
        step = rng/6
        ylim = (y_list[i].min() - step, y_list[i].max() + step)
        yticks = np.arange(y_list[i].min(), y_list[i].max()+0.0002, step=step*2)
        axs[i].set(ylim=ylim, yticks=yticks)

      # create x-axis ticks and labels
      axs[i].set(xlabel=xlabel)
      rng = (x.max() - x.min())
      step = rng/6
      xlim = (x.min() - step, x.max() + step)
      xticks = np.arange(x.min(), x.max()+.001, step=step*2)
      axs[i].set(xlim=xlim, xticks=xticks)

      # plot data
      for i in range(0, len(axs)):
        axs[i].scatter(x, y_list[i], color='thistle', s=30, edgecolors='k', linewidths=0.8,
             label='rms: ' + str(np.round(np.std(y_list[i]), 3)))

        leg = axs[i].legend(handlelength=0, handletextpad=0, loc='upper left')
        for item in leg.legendHandles:
          item.set_visible(False)


      # align y axis labels
      fig.align_ylabels(axs)

      # save figure
      fig_path = os.path.join(ImgDir.IMG_DIR, Inputs.fig_title)
      fig.savefig(fig_path, dpi=1200)


def process(files):
  good_files = []
  for file in files:
    name, extension = os.path.splitext(file)
    if extension == '.fits':
      good_files.append(file)
  map_seq=sunpy.map.Map(sorted(good_files))

  vmaps,mmaps,imaps=[],[],[]
  for j, map_obj in enumerate(map_seq):
    if map_obj.meta['content'] == 'DOPPLERGRAM':
      vmaps.append(map_obj)
    elif map_obj.meta['content'] == 'MAGNETOGRAM':
      mmaps.append(map_obj)
    elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
      imaps.append(map_obj)

  D={}
  for i in map_seq:
    t1=str(i.date)
    if t1 not in D:D[t1]=[i]
    else:D[t1].append(i)

  csvdata=[]
  if os.path.exists("/Users/suoenallecsim/Downloads/AllProcessedData.csv"):
    with open("/Users/suoenallecsim/Downloads/AllProcessedData.csv",'r') as file:
      L=file.read()
      L=[eval(i)for i in L.split('\n')[:-1]]
      for i in L:
        if i[0] in D:del D[i[0]]
  sys.stdout=open("/Users/suoenallecsim/Downloads/AllProcessedData.csv",'a')
  for date in D:
    maps=D[date]
    if len(maps)==0:pass
    x, y, pdim, r, d, mu = ctfuncs.coordinates(maps[0])
    wij, nij, rij = ctfuncs.vel_coords(x, y, pdim, r, maps[0])

    maps = ctfuncs.fix_mu(mu, maps, mu_cutoff=Parameters.mu_cutoff)
    vmap,mmap,imap=None,None,None
    for map_obj in maps:
      if map_obj.meta['content'] == 'DOPPLERGRAM':
        vmap = map_obj
      elif map_obj.meta['content'] == 'MAGNETOGRAM':
        mmap = map_obj
      elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
        imap = map_obj

    deltaw, deltan, deltar, dij = sfuncs.rel_positions(wij, nij, rij, maps[0])
    vsc = sfuncs.spacecraft_vel(deltaw, deltan, deltar, dij, maps[0])
    a_parameters = [Parameters.a1, Parameters.a2, Parameters.a3]
    vrot = sfuncs.solar_rot_vel(wij, nij, rij, deltaw, deltan, deltar, dij, maps[0], a_parameters)
    pixA_hem = ctfuncs.pix_area_hem(wij, nij, rij, maps[0])#takes long
    Lij = lbfuncs.limb_polynomial(maps[0])
    rv_model, v_quiet, v_disc, v_phot, v_conv, f_bright, f_spot, f, unsigned_obs_flux, vphot_bright,vphot_spot, f_small, f_large, f_network, f_plage, quiet_flux, ar_flux,conv_flux, pol_flux, pol_conv_flux, vconv_quiet, vconv_large, vconv_small, int_flux=[None]*24
    if vmap:
      corrected_vel = vmap.data - np.real(vsc) - np.real(vrot)
      map_vel_cor = sfuncs.corrected_map(corrected_vel, vmap, map_type='Corrected-Dopplergram',frame=frames.HeliographicCarrington)


    if imap:
      Iflat = imap.data / Lij
      map_int_cor = sfuncs.corrected_map(Iflat, imap, map_type='Corrected-Intensitygram',frame=frames.HeliographicCarrington)
      int_flux = np.nansum(imap.data)
    if mmap:
      Bobs, Br = sfuncs.mag_field(mu, mmap, B_noise=Parameters.B_noise, mu_cutoff=Parameters.mu_cutoff)
      map_mag_obs = sfuncs.corrected_map(Bobs, mmap, map_type='Corrected-Magnetogram',frame=frames.HeliographicCarrington)
      map_mag_cor = sfuncs.corrected_map(Br, mmap, map_type='Corrected-Magnetogram',frame=frames.HeliographicCarrington)
      active, quiet = sfuncs.mag_thresh(mu, mmap, Br_cutoff=Parameters.Br_cutoff, mu_cutoff=Parameters.mu_cutoff)
      area = sfuncs.area_calc(active, pixA_hem)

    if imap and mmap:
      fac_inds, spot_inds = sfuncs.int_thresh(map_int_cor, active, quiet)
      thresh_arr = sfuncs.thresh_map(fac_inds, spot_inds)
      map_full_thresh = sfuncs.corrected_map(thresh_arr, mmap, map_type='Threshold',frame=frames.HeliographicCarrington)
      f_bright, f_spot, f = sfuncs.filling_factor(mu, mmap, active, fac_inds, spot_inds, mu_cutoff=Parameters.mu_cutoff)

      f_small, f_large, f_network, f_plage = sfuncs.area_filling_factor(active, area, mu, mmap, fac_inds,athresh=Parameters.athresh,mu_cutoff=Parameters.mu_cutoff)

      v_phot, vphot_bright, vphot_spot = sfuncs.v_phot(quiet, active, Lij, vrot, imap, mu, fac_inds, spot_inds, mu_cutoff=Parameters.mu_cutoff)
      unsigned_obs_flux = sfuncs.unsigned_flux(map_mag_obs, imap)
      quiet_flux, ar_flux, conv_flux, pol_flux, pol_conv_flux = sfuncs.area_unsigned_flux(map_mag_obs, imap,area,active,athresh=Parameters.athresh)
    if vmap and imap:
      v_disc = sfuncs.v_disc(map_vel_cor, imap)
    if vmap and imap and mmap:
      v_quiet = sfuncs.v_quiet(map_vel_cor, imap, quiet)
      v_conv = v_disc - v_quiet
      vconv_quiet, vconv_large, vconv_small = sfuncs.area_vconv(map_vel_cor, imap, active, area, athresh=Parameters.athresh)
      rv_model = rvs.calc_model('NEID', v_conv, v_phot)

    save_vals = [str(date),rv_model, v_quiet, v_disc, v_phot, v_conv, f_bright, f_spot, f, unsigned_obs_flux, vphot_bright,vphot_spot, f_small, f_large, f_network, f_plage, quiet_flux, ar_flux,conv_flux, pol_flux, pol_conv_flux, vconv_quiet, vconv_large, vconv_small, int_flux]
    print(save_vals)
    csvdata.append(list(save_vals))
    #rounded = np.around(save_vals, 3)
    #round_vals = [date_str, date_jd]
    #for val in rounded:
    # round_vals.append(val)
    #utils.append_list_as_row(csv_file, round_vals)
  sys.stdout=normalstdout
  return csvdata

#import allProcessedData
#To combine csv's: "cat *.csv > ../combinedProcessedData.csv"
#no need to worry about duplicates or sorting
with open("/Users/suoenallecsim/Documents/combinedProcessedData.csv",'r') as file:
  L=file.read()
  L=L.replace('nan','None')
  L=[eval(i)for i in L.split('\n') if i.strip()]
  L=sorted(i for i in set(tuple(i)for i in L) if len(i)>2)# and type(i)==list)
  #L[0][0]='2021.01.11_12:57:2'

Ldates=[i[0]for i in L]
Ldatesjd=[]
for i in Ldates:
  try:t1=datetime.datetime.strptime(i[:10],"%Y.%m.%d")
  except:t1=datetime.datetime.strptime(i[:10],"%Y_%m_%d")
  Ldatesjd.append(t1.date().toordinal()+1721425+(int(i[11:13])+(int(i[14:16])+int(i[17])/6)/60)/24)
Ldatesjd.sort()
Lvdisc=[i[3]for i in L]
Lvquiet=[i[2]for i in L]
Lrvmodel=[i[1]for i in L]
Lfbright=[i[6]for i in L]
Lfspot=[i[7]for i in L]
Lvdiscx,Lvdiscy=list(zip(*[(Ldatesjd[i],Lvdisc[i])for i in range(len(Ldatesjd)) if Lvdisc[i]]))
Lvdiscx,Lvdiscy=list(Lvdiscx),list(Lvdiscy)
Lvquietx,Lvquiety=list(zip(*[(Ldatesjd[i],Lvquiet[i])for i in range(len(Ldatesjd)) if Lvquiet[i]]))
Lvquietx,Lvquiety=list(Lvquietx),list(Lvquiety)
Lrvmodelx,Lrvmodely=list(zip(*[(Ldatesjd[i],Lrvmodel[i])for i in range(len(Ldatesjd)) if Lrvmodel[i]]))
Lrvmodelx,Lrvmodely=list(Lrvmodelx),list(Lrvmodely)
Lfspotx,Lfspoty=list(zip(*[(Ldatesjd[i],Lfspot[i])for i in range(len(Ldatesjd)) if Lfspot[i]]))
Lfspotx,Lfspoty=list(Lfspotx),list(Lfspoty)
Lvdiscyma=movingAverage(20,Lvdiscy)


#ax1=plt.subplot(3,1,1);plt.plot(Lx,Ly);plt.subplot(3,1,2,sharex=ax1);plt.plot(Lvdiscx,Lvdiscy);plt.subplot(3,1,3,sharex=ax1);plt.plot(Lrvmodelx,Lrvmodely);plt.show()


Lrvmodel1=[];Lvdisc1=[];Lfspot1=[];Lfbright1=[];
t1=closestListValue(Ldatesjd[-1],Lx)
Lxtrunc=Lx[:t1]
Lytrunc=Ly[:t1]
Lymatrunc=Lyma[:t1]
for i in Lxtrunc:
    t1=closestListValue(i,Ldatesjd)
    
    Lrvmodel1.append(L[t1][1]/1000)
    Lvdisc1.append(L[t1][3])
    Lfspot1.append(L[t1][7])
    Lfbright1.append(L[t1][6])

Lx1trunc=[]
for i in Lx1:
  if i[0]>Lxtrunc[-1]:
    break
  Lx1trunc.append(i)
t1=0
for i in Lx1trunc[-1]:
  if i>Lxtrunc[-1]:
    t1=Lx1trunc[-1].index(i)
    break
Lx1trunc[-1]=Lx1trunc[-1][:t1]
Lxtrunci={Lxtrunc[i]:i for i in range(len(Lxtrunc))}
Lymatrunc1=[]
Lrvmodel11=[]
for i in Lx1trunc:
  t1=[]
  t3=[]
  for j in i:
    t1.append(Lymatrunc[Lxtrunci[j]])
    t3.append(Lrvmodel1[Lxtrunci[j]])
  t2=sum(t1)/len(t1)
  t4=sum(t3)/len(t3)
  for j in t1:
    Lymatrunc1.append(j-t2)
  for j in t3:
    Lrvmodel11.append(j-t4)

Lx1trunc0=[i[0]for i in Lx1trunc]
Dxytrunc={Lxtrunc[i]:Lytrunc[i]for i in range(len(Lxtrunc))}
Dxrvmodeltrunc={Lxtrunc[i]:Lrvmodel1[i]for i in range(len(Lxtrunc))}
Lymeans=[sum(Dxytrunc[j]for j in i)/len(i)for i in Lx1trunc]
Lrvmodelmeans=[sum(Dxrvmodeltrunc[j]for j in i)/len(i)for i in Lx1trunc]
Lymedians=[float(np.median([Dxytrunc[j]for j in i]))for i in Lx1trunc]
Lrvmodelmedians=[float(np.median([Dxrvmodeltrunc[j]for j in i]))for i in Lx1trunc]

#ax1=plt.subplot(2,1,1);plt.plot(Lxtrunc,Lytrunc);plt.subplot(2,1,2,sharex=ax1,sharey=ax1);plt.plot(Lxtrunc,Lrvmodel1);plt.show()
#plt.scatter(Lytrunc,Lrvmodel1);plt.show()
#t1=plt.hist2d(Lytrunc,Lrvmodel1,bins=(40,40));plt.show()
#plt.scatter(Lfspot1,Lytrunc);plt.show()
#plt.scatter(Lfspot1,Lrvmodel1);plt.show()
#plt.scatter(Lfbright1,Lytrunc);plt.show()
#plt.scatter(Lfbright1,Lrvmodel1);plt.show()
#scipy.stats.pearsonr(Lytrunc,Lfbright1);scipy.stats.pearsonr(Lrvmodel1,Lfbright1)
#ax1=plt.subplot(2,1,1);plt.plot(Lxtrunc,Lymatrunc1);plt.subplot(2,1,2,sharex=ax1);plt.plot(Lxtrunc,Lrvmodel11);plt.show()
#ax1=plt.subplot(2,1,1);plt.plot(Lx1trunc0,Lymeans);plt.subplot(2,1,2,sharex=ax1);plt.plot(Lx1trunc0,Lrvmodelmeans);plt.show()

#plt.plot(Ldatesjd[:7000],Lrvmodel[:7000]);plt.show()
#t1=np.fft.rfft(Lrvmodel[:7000]);plt.plot(range(len(t1)),np.abs(t1));plt.show()
#t1=np.fft.rfft(Lrvmodel);plt.plot(range(len(t1)),np.abs(t1));plt.show()


def removeOutliers(L,L2=[],n=1,bound=5):
    "n iterations"
    "L2 is an optional second list whose values are removed at the same indices as in L"
    "bound is a multiplier for the mean of differences to determine what values are removed"
    L=list(L);L2=list(L2)
    for i in range(n):
        t1=[0]+[abs(L[i]-L[i-1])for i in range(1,len(L))]
        t2=bound*np.mean(t1)
        L=[L[i]for i in range(len(t1))if t1[i]<t2]
        if L2:L2=[L2[i]for i in range(len(t1))if t1[i]<t2]
    if L2:return L,L2
    return L

def removeOutliers2(Ldatesjd0,Lrvmodel0,n=7):
  for k in range(n):
    ti=[]
    for i in range(len(Ldatesjd0)//2160-1):
      t1,t2=Ldatesjd0[2160*i:2160*(i+1)],Lrvmodel0[2160*i:2160*(i+1)];t3,t4=np.mean(t2),np.std(t2)
      for j in range(len(t2)):
        if not t3-4*t4<t2[j]<t3+4*t4:
          ti.append(2160*i+j)
    for i in range(len(ti)):
      del Ldatesjd0[ti[i]-i]
      del Lrvmodel0[ti[i]-i]
  return Ldatesjd0,Lrvmodel0


#remove outliers
Ldatesjd0=list(Ldatesjd);Lrvmodel0=list(Lrvmodel)

Ldatesjd0,Lrvmodel0=removeOutliers2(Ldatesjd0,Lrvmodel0,7)

Lrvmodel0,Ldatesjd0=removeOutliers(Lrvmodel0,Ldatesjd0,60)



#fill in missing data
def fillData(Ldatesjd0,Lrvmodel0,freq=[40,50]):
  t1=[]
  for i in range(1,len(Ldatesjd0)):
    if round((Ldatesjd0[i]-Ldatesjd0[i-1])*3600*24)not in freq:
      t1.append(Ldatesjd0[i-1]+int(np.mean(freq))/(3600*24))
      while round((Ldatesjd0[i]-t1[-1])*3600*24)>max(freq):
        t1.append(t1[-1]+int(np.mean(freq))/(3600*24))
  Ldatesjd0,Lrvmodel0=Ldatesjd0+t1,Lrvmodel0+[0 for i in t1]
  Ldatesjd0,Lrvmodel0=zip(*sorted(zip(Ldatesjd0,Lrvmodel0)))
  Ldatesjd0,Lrvmodel0=list(Ldatesjd0),list(Lrvmodel0)
  for i in range(1,len(Lrvmodel0)):
    if Lrvmodel0[i]==0:Lrvmodel0[i]=Lrvmodel0[i-1]
  return Ldatesjd0,Lrvmodel0

Ldatesjd0,Lrvmodel0=fillData(Ldatesjd0,Lrvmodel0)

#t1=savgol_filter(Lrvmodel0,192,1)#192 points is 1/10 of a day
Lrvmodel0sav=savgol_filter(Lrvmodel0,32,1)
Lrvmodel0savsav=savgol_filter(Lrvmodel0sav,32,1)
#plt.plot(Ldatesjd0,Lrvmodel0);plt.plot(Ldatesjd0,Lrvmodel0savsav);plt.show()#t1);plt.show()
Lrvmodelsub=[Lrvmodel0[i]-Lrvmodel0savsav[i]for i in range(len(Lrvmodel0))]
#plt.plot(Ldatesjd0,Lrvmodelsub);plt.show()

with open("/Users/suoenallecsim/Documents/rv_binned.csv",'r') as file:
    Lbin=[[float(j)for j in i.split(',') if j]for i in str(file.read()).split('\n')[1:]]
    Lbinx=[i[0]for i in Lbin if i]
    Lbiny=[i[1]for i in Lbin if i]
    #plt.plot(Lbinx,Lbiny);plt.show()

Lx0,Ly0=fillData(Lx,Ly,[82,83,93,94])
Ly0sav=savgol_filter(Ly0,16,1)
Ly0savsav=savgol_filter(Ly0sav,16,1)
#plt.plot(Lx0,Ly0);plt.plot(Lx0,Ly0savsav);plt.show()#t2);plt.show()
Lysub=[Ly0[i]-Ly0savsav[i]for i in range(len(Ly0))]
#plt.plot(Lx0,Lysub);plt.show()

