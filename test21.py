from test import *
normalstdout=sys.stdout
sys.stdout=open("multilog",'a')
print("TEST21 START SUCCESSFUL")
print(f"Current Time: {time.time()}\n")

#with open("current_process_num",'r') as file:
#    n=int(file.read())
n=sys.argv[1]
print(f"Process: {n}\n")
sys.stdout.close()
#print(sys.argv[1])

process(f"/storage/home/dcc5480/scratch/multi/{n}/jsoc/",f"/storage/home/dcc5480/work/Data/multi_{n}_{time.time()}.csv")

os.rename(f"/storage/home/dcc5480/scratch/multi/{n}/lastProcess",f"/storage/home/dcc5480/scratch/multi/{n}/lastProcess_done")

sys.stdout=open("multilog",'a')
print(f"TEST21 FINISHED! (process {n})")
print(f"Current Time: {time.time()}")



sys.stdout.close()


