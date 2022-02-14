import numpy as n
import re
import stuffr
import matplotlib.pyplot as plt
import scipy.constants as sc
import datetime as dt
import scipy.interpolate as sio
monthnums={"Feb":2,"Jan":1,"Mar":3}

class horizons_moon:
    def __init__(self,fname="horizons_2022_02_13.txt",ipp=32e-3):
        self.Rmoon=1737.1e3
        f=open(fname,"r")
        soe=False
        self.r=[]
        self.ipp=ipp
        self.dts=[]        
        self.t=[]
        self.rr=[]
        for l in f.readlines():
            if l.strip() == "$$EOE":
                soe=False
            if soe:
                toks=l.strip().split(",")
                datestr=toks[0]
                range_km = float(toks[11])
                range_rate_kms=float(toks[12])
                reg=re.search("(....)-(...)-(..) (..):(..):(..)",datestr)
                year=int(reg.group(1))
                month=monthnums[reg.group(2)]
                day=int(reg.group(3))
                hour=int(reg.group(4))
                mins=int(reg.group(5))
                sec=int(reg.group(6))
                ut=stuffr.date2unix(year,month,day,hour,mins,sec)
                self.dts.append(dt.datetime.fromtimestamp(ut))
                self.r.append(range_km)
                self.rr.append(range_rate_kms)
                self.t.append(ut)
                
#                print("%s %f %f"%(stuffr.unix2datestr(ut),range_km,range_rate_kms))
            
            if l.strip() == "$$SOE":
                soe=True
        self.r=n.array(self.r)
        self.rr=n.array(self.rr)
        self.t=n.array(self.t)
        self.rfun=sio.interp1d(self.t,self.r)
        self.rrfun=sio.interp1d(self.t,self.rr)
        hdt=self.t[1]-self.t[0]
        self.rr2=n.gradient(self.r,hdt)
        

        
        
    def plot_range(self):
        plt.plot(self.t,self.r)
        plt.show()

        # leading edge
        plt.plot(self.dts,1e3*n.mod(2.0*(self.r*1e3-self.Rmoon)/sc.c,self.ipp))
        # limb
        plt.plot(self.dts,1e3*n.mod(2.0*(self.r*1e3)/sc.c + 1.69e-3,self.ipp))        
        plt.axhline(1.69)
        plt.show()
        


if __name__ == "__main__":
    hm=horizons_moon()
    hm.plot_range()
