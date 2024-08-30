# NEID-Research
Absolute filepaths are used for many cases, so they will need to be modified. For downloading and processing data, the "scratch" directory of the ROAR Collab supercomputer is used, as each 3-day group of data takes around 300 GB of space.

## NEID-Research.py
Used for experimenting with data from the NEID telescope and processed STO data. See comments in file for details.

## test.py
Includes functions:
- fullPipeline (Removed, see older commits: old, not used anymore, see process1)
  - Converts downloaded .fits files into data including modeled rv's. Is able to convert times without all three fits files (magnetogram, dopplergram, intensitygram) into partial data output.
- fidoSearch
  - Searches for and downloads .fits files. Retries download where errors occur (won't retry successful downloads).
- process
  - Checks for times with enough data and sends the .fits files to process1.
- process1
  - Runs the pipeline to convert .fits files into rv model data. Prints processed data to file with variable name "outpath", deletes used .fits files, and prints errors to ProcessingErrors.csv.
    - outpath is determined by the julian date (first day of a 3-day group), passed from multiProcess to main_job.slurm, then to test21.py where outpath is set to a csv file with name of the form "multi_{day}_{random-int}.csv", then to process, and to process1 where the data is written.
- multiMultiProcess
  - Finishes incomplete processes, using "lastProcess" files from multiProcess and writes the current time to them again in case the process doesn't finish again.
  - Checks file "multiMultiNums" for what day to start from. By default, the day to process to is 10 days before the present. Writes to multiMultiNums after reading, by default 15 days later than what was there before (determined by variable "firstnum").
  - Starts new multiProcess function (for 15 days of data by default, also determined by "firstnum").
  - Runs multi_job.slurm which will start a new job to run multiMultiProcess again.
  - Each job can only run up to 48 hours, and multiProcess can be run (at least) ~5 times per job, so a new job is requested for every five runs of multiProcess.
- multiProcess
  - Inputs a date range and splits it into groups of 3 days, rounding up. This amount of days is used so that .tar files are downloaded, which is faster than individual .fits files.
  - Attempts to find and download STO data for each set of three days, will try 4+ days if there isn't enough data.
  - Unpacks and deletes .tar files
  - Writes the time to file "lastProcess" to indicate when it started processing the .fits files data. Each such file is for a group of 3 days, and it takes longer than 48 hours to process one of these days (more than 1 job), so each time multiMultiProcess is run, it checks if each "lastProcess" file has a time from more than 48 hours ago, and if so tries to continue processing its respective 3-day group.
  - Runs main_job.slurm to process each 3-day group.

## multi_job.slurm
Runs multi.py in a new job, is run from multiMultiProcess after it finishes running multiProcess. Only one of these jobs can run at a time.

## multi.py
Runs multiMultiProcess in a new job, is run from multiProcess after it downloads and unpacks data.

## main_job.slurm
Runs test21.py in a new job. Many of these jobs can run in parallel.

## test21.py
Runs process on a 3-day group. If it finishes running process, then it renames its respective "lastProcess" file to "lastProcess_done" to indicate that multiMultiProcess doesn't need to try to continue processing this group. 

## log (example file included here)
Includes info from fidoSearch, including days downloaded and errors. This file tends to be large.

## multilog (example file included here)
Includes timestamps for when test21.py, multi.py, multiMultiProcess, and multiProcess start and finish. Generally best for troubleshooting issues.
