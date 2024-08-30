# NEID-Research

## NEID-Research.py
Used for experimenting with data from the NEID telescope and processed STO data.

## test.py
Includes functions:
- fullPipeline (Old, not used anymore)
  - Converts downloaded .fits files into data including modeled rv's.
- fidoSearch
  - Searches for and downloads .fits files.
- process
  - Checks for times with enough data and sends the .fits files to process1.
- process1
  - Runs the pipeline to convert .fits files into rv model data.
- multiMultiProcess
  - Finishes incomplete processes, using "lastProcess" files from multiProcess and writes the current time to them again in case the process doesn't finish again.
  - Starts new multiProcess functions.
  - Runs multi_job.slurm which will start a new job to run multiMultiProcess again.
  - Each job can only run up to 48 hours, and multiProcess can be run (at least) ~5 times per job, so a new job is requested for every five runs of multiProcess.
- multiProcess
  - Inputs a date range and splits it into groups of 3 days. This amount of days is used so that .tar files are downloaded, which is faster than individual .fits files.
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


