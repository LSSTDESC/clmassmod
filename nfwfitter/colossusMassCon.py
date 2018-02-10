#########
# Defines mass-con relations
#########

import numpy as np
import nfwutils
import nfwmodeltools as tools
import colossus.cosmology.cosmology as cCosmo
import colossus.halo.concentration as chc


class IterateC200(object):

    def __call__(self, m, z, overdensity = 200, tolerance = 0.05):

        """Iterate to calculate c200 from an M_delta that isn't m200.

        m           - in M_sun/h
        z           - redshift
        overdensity - wrt to the critical density

        Inspired by J. Dietrich, 2013
        Heavily optimized by D. Applegate
        """

        rho_crit = nfwutils.global_cosmology.rho_crit(z)

        c0 = self.mc(m, z)

        if overdensity == 200:
            return c0

        oldC = 100.
        while (np.abs(c0 - oldC)/oldC >= 0.05):
            oldC = c0
            rscale = tools.rscaleConstM(m, c0, rho_crit, overdensity)
            r200 = c0*rscale
            m200c = 200*rho_crit*(4*np.pi/3)*r200**3
            c0 = self.mc(m200c,z)

        return c0


#########

def matchCosmo():

    curcosmo = nfwutils.global_cosmology
    isMatch = False
    try:
        cur_ccosmo = cCosmo.getCurrent()

        isMatch = curcosmo.H0 == cur_ccosmo.H0 and \
                  curcosmo.omega_m == cur_ccosmo.Om0 and \
                  curcosmo.omega_l == cur_ccosmo.OL0
    except:
        pass


    if isMatch is False:

        cCosmo.setCosmology('curCosmo', dict(H0 = curcosmo.H0,
                                             Om0 = curcosmo.omega_m,
                                             Ode0 = curcosmo.omega_l,
                                             Ob0 = 0.049,
                                             sigma8 = 0.9,
                                             ns = 0.95))



#########


class ColossusMC(IterateC200):

    def configure(self, config):

        self.modelname = config['colossusmcname']


    ###


        


    ###


    def mc(self, m, z):

        matchCosmo()

        c200 = chc.concentration(m, '200c', z, self.modelname)

        return c200


        

        


matchCosmo()

