from test import *
normalstdout=sys.stdout
sys.stdout=open("multilog",'a')
print("TEST2 START SUCCESSFUL")
print(f"Current Time: {time.time()}\n")

with open("current_process_num",'r') as file:
    n=int(file.read())
print(f"Process: {n}\n")
sys.stdout.close()

process(f"/storage/home/dcc5480/scratch/multi/{n}/jsoc/",f"/storage/home/dcc5480/work/Data/multi_{n}_{time.time()}.csv")

sys.stdout=open("multilog",'a')
print(f"TEST2 FINISHED! (process {n})")
print(f"Current Time: {time.time()}")



sys.stdout.close()


