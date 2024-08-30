from test import *
normalstdout=sys.stdout
sys.stdout=open("multilog",'a')
print("MULTI START SUCCESSFUL")
print(f"Current Time: {time.time()}\n")

with open("multiMultiNums",'r') as file:
    n1,n2=file.read().split(' ')
    n1,n2=int(n1),int(n2)
print(f"multiMultiProcess: {n1} {n2}\n")
sys.stdout.close()

multiMultiProcess(n1,n2)
#process(f"/storage/home/dcc5480/scratch/multi/{n}/jsoc/",f"/storage/home/dcc5480/work/Data/multi_{n}_{time.time()}.csv")

sys.stdout=open("multilog",'a')
print(f"MULTI FINISHED! (process {n1} {n2})")
print(f"Current Time: {time.time()}")



sys.stdout.close()


