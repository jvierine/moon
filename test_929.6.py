import glob
import numpy as n
import scipy.io as sio
import matplotlib.pyplot as plt
import scipy.constants as c
import os

if __name__ == "__main__":

    # sshfs eiscat@goppi.eiscat.uit.no:/ goppi
    fl=glob.glob("/home/j/goppi/data/dimmer1_1.0_NO@uhf/*/*.mat")
    fl.sort()

    # make sample index
    idx=n.arange(32001,dtype=n.float64)

    # store phase here

    # continuous ipp
    zipp=n.zeros(32000,dtype=n.complex64)
    pha0=0.0
    # look at last 10 files
    nback=10
    for i in range(nback):
        fname=fl[len(fl)-nback-1+i]
        d=sio.loadmat(fname)
        z=n.array(d["d_raw"][:,0],dtype=n.complex64)

        zmed=[]
        # this many ipps in file
        n_ipp=int(len(z)/31930)
        for k in range(n_ipp):
            print(k)
            zipp[60:(60+31930)]=z[ (k*31930) : (k*31930 + 31930) ]
#            csin=n.exp(-1j*2.0*n.pi*400e3*idx/1e6)*n.exp(-1j*pha0)
            # angle 
 #           pha0=n.angle(csin[len(csin)-1])
            zd=zipp#*csin[0:len(zipp)]
            zmed.append(n.median(zd))
        zmed=n.array(zmed)
        plt.title(fname)
        plt.plot(zmed.real)
        plt.plot(zmed.imag)
        plt.show()

