import glob
import numpy as n
import scipy.io as sio
import matplotlib.pyplot as plt
import stuffr
import scipy.constants as c
import digital_rf as drf
import os

class e2drf:
    def __init__(self,
                 dirname="/home/j/goppi/data/dimmer1_1.0_NO@uhf/*/*.mat",
                 sr=1000000,
                 outdirname="/data0/eiscat/moon/2022_02_14/uhf",
                 ipp=32000,
                 L=31930,
                 drawlen=6686000,
                 pad=60):

        fl=glob.glob(dirname)
        fl.sort()

        # look for files that are good 
        #fl2=[]
        #for f in fl:
        #    d=sio.loadmat(f)
        #    if len(d["d_raw"]) == drawlen:
        #        print(f)
        #        fl2.append(f)
        #    else:
        #        print("bad %s"%(f))
        #fl=fl2
        self.drawlen=drawlen
        self.fl=fl
        self.L=L
        self.pad=pad
        self.ipp=ipp
        self.sr=sr

        os.system("mkdir -p %s"%(outdirname))

        d=sio.loadmat(fl[0])
        if len(d["d_raw"]) != drawlen:
            print("wrong vector len!")
            exit(0)
        
        # samples since 1970
        t0=int(self.get_t0(d))*self.sr
        print(t0)
        
        self.w=drf.DigitalRFWriter(outdirname, n.complex64, 3600, 1000, t0, 1000000, 1, "32m", compression_level=0, checksum=False, is_complex=True, num_subchannels=1, is_continuous=True, marching_periods=True)


    def get_t0(self,d):
        # we need to work a bit to get the exact start time.
        dpar=d["d_parbl"][0]
        dt=dpar[5]-int(dpar[5])
        # this is the approxiamte time of the end of the dump
        dump_end=stuffr.date2unix(int(dpar[0]),int(dpar[1]),int(dpar[2]),int(dpar[3]),int(dpar[4]),int(dpar[5]))+dt
        # this is the length of the dump
        dump_dt=dpar[6]
        self.dump_dt=dump_dt
        # this many files since experiment start
        didx=int(dpar[11]-1)
        # this is the exact start time of the experiment
        exp_start=n.round(dump_end - didx*dump_dt)
        
        # exact start time of this file
        t0 = exp_start + didx*dump_dt
        return(t0)
        
    def map_files(self):
        zipp = n.zeros(self.ipp,dtype=n.complex64)
        L=self.L

        for f in self.fl:
            d=sio.loadmat(f)
            z=n.array(d["d_raw"],dtype=n.complex64)
            if len(z) != self.drawlen:
                raise Exception("wrong vector length!")
            print(stuffr.unix2datestr(self.get_t0(d)))
            n_ipp=int(self.dump_dt/(self.ipp/1e6))

            for i in range(n_ipp):
                zipp[self.pad:(self.pad+L)]=z[ (L*i):(L*i+L), 0 ]
                self.w.rf_write(zipp)
                


if __name__ == "__main__":

#     
#    m=e2drf(dirname="/home/j/goppi/data/dimmer1_1.0_NO@uhf/",
#            outdirname="/data0/eiscat/moon_2022_02_13")                  
#    m.map_files()

    m=e2drf(dirname="/home/j/goppi/data/dimmer2_1.0_NO@uhf/20220214_17/*.mat",
            outdirname="/data0/eiscat/moon/2022.02.14/uhf",
            ipp=33500,
            L=33430,
            drawlen=6686000,
            pad=60)

    m.map_files()
    
    
    
