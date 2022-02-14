import glob
import numpy as n
import scipy.io as sio
import matplotlib.pyplot as plt
import horizons as hor
import stuffr
import scipy.constants as c
import digital_rf as drf

class moon_mapper:
    def __init__(self,
                 eph,
                 dirname="/data0/eiscat/",
                 ch="moon_2022_02_13",
                 pad=3000,
                 L=31930,
                 ipp=32000,
                 bitlen=2,
                 txlen=1800,
                 N_ipp=200,
                 dec=10):
        
        d=drf.DigitalRFReader(dirname)
        b=d.get_bounds(ch)
        self.b=b
        self.d=d
        self.ch=ch
        self.eph=eph
        self.bitlen=bitlen
        self.L=L
        self.pad=pad
        self.ipp=ipp
        self.txlen=txlen
        self.dec=dec
        self.sr=1000000.0
        self.n_rg = int(1e6*2.0*self.eph.Rmoon/c.c + 2*pad)
        self.E = n.zeros([N_ipp,self.n_rg],dtype=n.complex64)
        self.T = n.zeros([N_ipp,self.n_rg],dtype=n.complex64)        
        self.S = n.zeros([N_ipp,self.n_rg],dtype=n.float32)
        self.nrd=int(self.n_rg/dec)
        self.SD = n.zeros([N_ipp,self.nrd],dtype=n.float32)


        self.N_ipp=N_ipp
        self.N_maps = int(n.floor( (self.b[1]-self.b[0])/(ipp*N_ipp)) )
        print(self.N_maps)


    def map_files(self):

        zipp = n.zeros(self.ipp,dtype=n.complex64)
        ztx = n.zeros(self.n_rg,dtype=n.complex64)        
        tv=n.array(n.arange(self.ipp+1)/1e6,dtype=n.float32)

        i0 = self.b[0]
        ipp_idx = 0
        ph0=0.0
        
        for mi in range(self.N_maps):
            img_t0=(self.b[0]+mi*self.N_ipp*self.ipp)/self.sr

            for ip in range(self.N_ipp):
                i0 = self.b[0] + ipp_idx*self.ipp
                ipp_idx+=1                
                t0=i0/self.sr

                rnow = 1e3*self.eph.rfun(t0) - self.eph.Rmoon
                rtt_s = 2.0*rnow/c.c
                rttnow = int(self.sr*rtt_s)
                
                ztx[:]=self.d.read_vector_c81d(i0,self.n_rg,self.ch)
                ztx[self.txlen:self.n_rg]=0.0

                print(t0 + rtt_s/2.0)
                
                rr=1e3*self.eph.rrfun(t0 + rtt_s/2.0)
                dopfreq=2.0*929.6e6*rr/c.c

                print("rtt %d dop %1.2f"%(rttnow,rr))
                csin=n.exp(-1j*2.0*n.pi*dopfreq*tv)*n.exp(1j*ph0)
                # save phase
                ph0=n.angle(csin[len(csin)-1])


#                ztx=ztx*n.exp(-1j*n.angle(ztx[88]))                
                self.T[ip,:]=ztx

                
                z_tx2 = ztx*csin[0:len(ztx)]

#                z_tx2 = csin[0:len(ztx)]                

                z_echo=self.d.read_vector_c81d(i0+rttnow-self.pad,self.n_rg,self.ch)

                zd=n.fft.ifft(n.fft.fft(z_echo)*n.conj(n.fft.fft(z_tx2)))
                
                self.E[ip,:]=zd


            if False:
                
                plt.plot(self.T[:,88].real)
                plt.plot(self.T[:,88].imag)
                plt.show()
                plt.subplot(121)
                plt.title("Transmit pulse phase")
                plt.pcolormesh(n.angle(self.T[:,0:1000]),cmap="hsv")
                plt.colorbar()
                plt.xlabel("Time (samples)")
                plt.ylabel("IPP")                
                #            plt.show()
                plt.subplot(122)
                plt.title("Echo phase")                
                plt.pcolormesh(n.angle(self.E[:,2000:4000]),cmap="hsv")
                plt.xlabel("Time (samples)")
                plt.ylabel("IPP")                
                
                plt.colorbar()
                plt.show()
            
  
            for ri in range(self.n_rg):
                self.S[:,ri]=n.abs(n.fft.fftshift(n.fft.fft(self.E[:,ri])))**2.0

            for fi in range(self.N_ipp):
                self.SD[fi,:]=stuffr.decimate(self.S[fi,:],dec=self.dec)[0:self.nrd]

            dB=10.0*n.log10(n.transpose(self.SD))
            nfloor=n.nanmedian(dB)
            dB=dB-nfloor
            print(mi)
            fvec=n.fft.fftshift(n.fft.fftfreq(self.N_ipp,d=self.ipp/1e6))
            rvec=(n.arange(self.nrd)-self.pad/self.dec)*0.15*self.dec
            plt.pcolormesh(fvec,rvec,dB,cmap="gray",vmin=-3,vmax=50.0)
            plt.title(stuffr.unix2datestr(img_t0))
            plt.xlabel("Doppler shift (Hz)")
            plt.ylabel("Sub-radar point range offset (km)")
            cb=plt.colorbar()
            cb.set_label("SNR (dB)")
            plt.tight_layout()
#            plt.show()            
            plt.savefig("plots/moon-%d.png"%(int(img_t0)))
            plt.clf()
            plt.close()

                
        








if __name__ == "__main__":
    if False:
        # 2022_02_13
        heph=hor.horizons_moon(fname="horizons_2022_02_13.txt",
                               dirname="/data0/eiscat/",
                               ch="moon_2022_02_13",
                               )
        m=moon_mapper(heph)
        m.map_files()

    if True:
        # 2022_02_14
        heph=hor.horizons_moon(fname="horizons_2022_02_13.txt")
        m=moon_mapper(heph,
                      dirname="/data0/eiscat/moon/2022.02.14",
                      ch="uhf",
                      ipp=33500,
                      L=33430)
        m.map_files()


