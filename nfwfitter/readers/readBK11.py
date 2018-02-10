######################
# Reads in BK11 files so that nfwfit.py can use them
######################

import astropy.io.fits as pyfits, ldac
import nfwutils
import numpy as np
import catalog

class BK11SimReader(object):

    def __init__(self, *args, **keywords):

        pass


    def getCosmology(self):
        return nfwutils.Cosmology(omega_m = 0.27, omega_l = 0.73, h = 0.7)

    def load(self, filebase):

        return BK11Sim(filebase)


class BK11Sim(catalog.Catalog):

    def __init__(self, filename):

        super(BK11Sim, self).__init__()


        sim = ldac.LDACCat(pyfits.open(filename)[1])

        boxlength = sim['BOXWIDTHCOMOVINGHINVMPC'] / nfwutils.global_cosmology.h
        gridsize = np.sqrt(float(sim['A00'].shape[1]))

        clusterz = float(sim['ZLENS'])
        self.zcluster = clusterz

        ####becker angles

        Ng = 512

        dc_readout = nfwutils.global_cosmology.comovingdist(clusterz)+sim['BOXLENGTHCOMOVINGHINVMPC'][0]/(2.0*nfwutils.global_cosmology.h)
        da = np.arctan2(sim['BOXWIDTHCOMOVINGHINVMPC'][0]/nfwutils.global_cosmology.h,dc_readout)/np.pi*180.0*60.0/Ng

        xcen = Ng/2.0
        ycen = Ng/2.0

        radii_arcmin = np.zeros((Ng, Ng))
        cosphi = np.zeros((Ng,Ng))
        sinphi = np.zeros((Ng,Ng))
        x_arcmin = np.zeros((Ng,Ng))
        y_arcmin = np.zeros((Ng,Ng))

        for i in range(Ng):
            for j in range(Ng):
                xr = (i + 0.5 - xcen)*da
                yr = (j + 0.5 - ycen)*da
                r = np.sqrt(xr*xr + yr*yr)
                radii_arcmin[i,j] = r
                x_arcmin[i,j] = xr
                y_arcmin[i,j] = yr

                cosphi[i,j] = xr/r
                sinphi[i,j] = yr/r


        A00 = sim['A00'][0]
        A11 = sim['A11'][0]
        A10 = sim['A10'][0]
        A01 = sim['A01'][0]

        kappa = 0.5*(2- A00 - A11)
        gamma1 = 0.5*(A11 - A00)
        gamma2 = -0.5*(A01+A10)


        self.x_arcmin = x_arcmin.flatten()
        self.y_arcmin = y_arcmin.flatten()

        Dl = nfwutils.global_cosmology.angulardist(clusterz)
        self.x_mpc = (self.x_arcmin/60.)*(np.pi/180.)*Dl
        self.y_mpc = (self.y_arcmin/60.)*(np.pi/180.)*Dl


        redshifts = np.ones(len(self.x_arcmin))*float(sim['ZSOURCE'])
        beta_s = nfwutils.global_cosmology.beta_s(redshifts, clusterz)

        self.gamma1_inf = gamma1.flatten()/beta_s
        self.gamma2_inf = gamma2.flatten()/beta_s
        self.kappa_inf  = kappa.flatten()/beta_s


        self.m500 = sim['M500C'][0]
        self.m200 = sim['M200C'][0]
        self.c200 = sim['C200C'][0]


