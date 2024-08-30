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
normalstdout=sys.stdout
sys.path.append(os.path.realpath('../../'))
import SolAster.tools.rvs as rvs
import SolAster.tools.utilities as utils
from SolAster.tools.settings import *
from SolAster.tools.plotting_funcs import hmi_plot
import traceback
import threading


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

    csv_name=Inputs.csv_name#start_date, end_date, cadence, csv_name = utils.check_inputs(CsvDir.CALC,  start_date,  end_date, cadence,  Inputs.csv_name)


    # print out csv title
    print("Beginning calculation of values for csv file: " + csv_name)

    # List of header strings
    row_contents = ['date_obs', 'date_jd', 'rv_model', 'v_quiet', 'v_disc', 'v_phot', 'v_conv', 'f_bright', 'f_spot', 'f',
        'Bobs', 'vphot_bright', 'vphot_spot', 'f_small', 'f_large', 'f_network', 'f_plage',
        'quiet_flux', 'ar_flux', 'conv_flux', 'pol_flux', 'pol_conv_flux', 'vconv_quiet', 'vconv_large',
        'vconv_small']

    # create file names
#    csv_file = os.path.join(CsvDir.CALC, csv_name)
#    bad_dates_csv = os.path.join(CsvDir.CALC, csv_name[:-4]+'_bad_dates.csv')
#    print(bad_dates_csv)
#    utils.append_list_as_row(csv_file, row_contents)


    # get hmi data products
    time_range = datetime.timedelta(seconds=22)
    physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

    # get dates list
    xy = (end_date - start_date).seconds + (end_date - start_date).days * 24 * 3600
    dates_list = [start_date + datetime.timedelta(seconds=cadence*x) for x in range(0, int(xy/cadence))]
    #remove the equal sign here to use this part
    for i, date in enumerate(dates_list):
    # convert the date to a string -- required for use in csv file
        date_obj, date_jd =date,date.toordinal()+1721425
        print(date_obj,date_jd)

        # pull image within specified time range
        result = Fido.search(a.Time(str(date_obj - time_range), str(date_obj + time_range)),
                         a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2])

        # add file to list
        file_download = Fido.fetch(result)
    return file_download
"""
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
            print('Calculations and save to file complete for ' + date_str + ' index: ' + str(i))"")
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

"""
def fidoSearch(start_date,end_date,path=None,max_conn=5):
    start_date=datetime.datetime(1858,11,17)+timedelta(start_date-2400000.5)
    end_date=datetime.datetime(1858,11,17)+timedelta(end_date-2400000.5)
    result=Fido.search(a.Time(str(start_date),str(end_date)),a.Instrument.hmi,a.Physobs.los_velocity | a.Physobs.los_magnetic_field | a.Physobs.intensity)
    sys.stdout=open("log",'a')
    print(time.time())
    print(result)
    #https://groups.google.com/g/sunpy/c/T9zVzvh0Iro?pli=1
    result1=Fido.fetch(result,path=path,max_conn=max_conn,overwrite=False)
    while len(result1.errors)>0:
        print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
        print(time.time())
        print(result1)
        print(result1.errors)
        result1=Fido.fetch(result1,path=path,max_conn=5,overwrite=True)
    sys.stdout.close()

def multiMultiProcess(start_date,end_date=int(time.time()/(86400)+2440577),firstnum=5):
    sys.stdout=open("multilog",'a')
    print(f"\nmultiMultiProcess Started: {start_date} {end_date} {time.time()}")
    print("\nContinuing Unfinished Processes...")
    s="/storage/home/dcc5480/scratch/multi/"
    for i in os.listdir(s):
        print(f"\nChecking {i}... {time.time()}")
        s1=s+i+"/jsoc/"
        if "lastProcess" not in os.listdir(s+i):
            if "lastProcess_done" not in os.listdir(s+i):print(f"Error, no lastProcess file: {i}")
        elif not os.path.exists(s1):print(f"Error, no jsoc folder: {i}")
        elif len(os.listdir(s1))>5:
            s1=s+i+"/lastProcess"
            with open(s1,'r') as file:t1=float(file.read())
            if time.time()-t1>173000:
                print(f"\nProcessing {i} {time.time()}")
                with open(s1,'w') as file:file.write(str(time.time()))
                os.system(f"sbatch --export=TEST_DATE='{i}' main_job.slurm")
    print("Preprocessing Done")
    print(f"Starting multiProcess {time.time()}")

    sys.stdout.close()
    with open("/storage/home/dcc5480/work/multiMultiNums",'w') as file:
        file.write(f"{start_date+firstnum*3} {end_date}")
    multiProcess(start_date,start_date+firstnum*3-3)
    os.system("sbatch multi_job.slurm")
    sys.stdout=open("multilog",'a')
    print(f"multi_job.slurm started {start_date} {end_date} {time.time()}")
    sys.stdout.close()


def multiProcess(start_date,end_date):
    sys.stdout=open("multilog",'a')
    print("\n"*10)
    print(start_date,end_date,time.time())
    for i in range(int(start_date),int(end_date)+1,3):#was +10#then +1,+3 allows tars
        sys.stdout.close()
        try:fidoSearch(i,i+3,f"/storage/home/dcc5480/scratch/multi/{i}/")
        except:
            sys.stdout=open("multilog",'a')
            print("\n\nDOWNLOAD QUANTITY ERROR: retrying fidoSearch with an extra day!")
            print(f"DAY: {i}\n\n")
            sys.stdout.close()
            fidoSearch(i,i+4,f"/storage/home/dcc5480/scratch/multi/{i}/")
            sys.stdout=open("multilog",'a')
            print("\n\nQUANTITY ERROR RESOLVED!")
            print(f"DAY: {i}\n\n")
            sys.stdout.close()
        os.system(f'cd /storage/home/dcc5480/scratch/multi/{i};for file in /storage/home/dcc5480/scratch/multi/{i}/*.tar; do tar -xf "$file";rm "$file"; done')
        #thread=threading.Thread(target=process,args=("/storage/home/dcc5480/scratch/multi/{i}/jsoc/",f"multi_{i}.csv"))#removed "jsoc/" since 1 day at a time instead of 10 does individual files instead of .tar#added again for 3
        #thread.start()
        #with open("current_process_num",'w') as file:
        #    file.write(str(i))
        os.system(f"sbatch --export=TEST_DATE='{i}' main_job.slurm")
        with open(f"/storage/home/dcc5480/scratch/multi/{i}/lastProcess",'w') as file:
            file.write(str(time.time()))
        sys.stdout=normalstdout
        print(f"Process {i} started...")
        sys.stdout=open("multilog",'a')
        print("Process Started",i,time.time())
    print("DONE")
    sys.stdout.close()

def process(files,outpath=None,skipfinished=True,delete=True):
        if type(files)==str:
                if files[-1]!='/':files+='/'
                files=[files+i for i in os.listdir(files)]
        if outpath==None:outpath="AllProcessedData.csv"
        good_files = []
        for file in files:
                name, extension = os.path.splitext(file)
                if extension == '.fits':
                        good_files.append(file)

        D={}
        #print(good_files)
        for i in good_files:
                t1=i[i.index('/')+1:]
                t1=t1[t1.index('45s')+4:]
                t1=t1[:t1.lower().index('_tai')-1]
                if t1 not in D:D[t1]=[i]
                else:D[t1].append(i)

        csvdata=[]
        if os.path.exists(outpath):
                with open(outpath,'r') as file:
                        L=file.read()
                        L=[eval(i)for i in L.split('\n') if i and 'nan' not in i]
                        for i in L:
                                if i[0] in D and skipfinished and not (None in i and len(D[i[0]])==3):del D[i[0]]
        process1(D,outpath,delete)

def process1(D,outpath,delete):
        csvdata=[]
        for date in D:
            try:
                    if len(D[date])<3:continue
                    sys.stdout=open(outpath,'a')
                    maps=D[date]
                    if len(maps)==0:pass
                    try:maps=sunpy.map.Map(sorted(maps))
                    except:
                            print([date])
                            pass
                    if len(D[date])==1:maps=[maps]
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
                    #	round_vals.append(val)
                    #utils.append_list_as_row(csv_file, round_vals)
                    sys.stdout.close()
                    if delete:
                        for i in D[date]:os.unlink(i)
            except Exception:
                sys.stdout.close()
                sys.stdout=open("ProcessingErrors.csv",'a')
                print([traceback.format_exc(),date,time.time()])
                sys.stdout.close()
            sys.stdout=normalstdout
        #return csvdata
