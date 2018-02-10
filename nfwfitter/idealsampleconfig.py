import readers.readMXXL
simreader = readers.readMXXL.MXXLSimReader()

import nfwfit
model = nfwfit.NFW_MC_Model()
massprior='linear'

import colossusMassCon
massconRelation = colossusMassCon.ColossusMC()
colossusmcname='diemer15'

import nfwfit
fitter = nfwfit.PDFScanner()

import rescalecluster
rescalecluster = rescalecluster.RedshiftRescaler()
targetz=0.9

import galaxypicker
galaxypicker = galaxypicker.AllGalaxyPicker()


import betacalcer
betacalcer = betacalcer.FixedBeta()
beta = 0.3


import shearnoiser
shearnoiser = shearnoiser.NoNoise()

import centergenerator
centergenerator = centergenerator.NoOffset()


import basicBinning
binner = basicBinning.BootstrapFixedBins()
profileMax = 1.5
profileMin = 0.5
binspacing = 'linear'
nbins = 12
profilecol = 'r_mpc'


import binnoiser
binnoiser = binnoiser.NoBinNoise()

import profilebuilder
profilebuilder = profilebuilder.ProfileBuilder()

