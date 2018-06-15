# set MaxEvents 10000000
#set RandomSeed 123
#
#  Phase II - Pile-Up
#
#  Main authors: Michele Selvaggi (UCL)
#
#  Released on:
#
#  Version: v02 beta - test TrackSmearing and vertexing
#
#
#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  RHadronFilter
  RHadronMerger

  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency
  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger
  TrackSmearing

  TimeSmearing

  VertexFinder

  TreeWriter
}

#################################
# Get stable R-Hadrons
#################################

module PdgCodeFilter RHadronFilter {
  set InputArray Delphes/allParticles
  set OutputArray stableRHadrons
  set Invert false

  set RequireStatus true
  set Status 104
}

#################################
# Merge R-Hadrons with other sable
#################################

module Merger RHadronMerger {
  add InputArray RHadronFilter/stableRHadrons
  add InputArray Delphes/stableParticles
  set OutputArray stableParticles
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray RHadronMerger/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.0

  # magnetic field
  set Bz 3.8
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  ## particles after propagation
  set InputArray  ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
	  (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
	  (abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.87) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.82) + \
	  (abs(eta) > 4.0) * (0.00)
  }
}


#####################################
# Electron tracking efficiency - ID
####################################

module Efficiency ElectronTrackingEfficiency {
  set InputArray  ParticlePropagator/electrons
  set OutputArray electrons
  # tracking efficiency formula for electrons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
	  (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
	  (abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0 && pt <= 10.0) * (0.82+pt*0.01) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 10.0) * (0.90) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt <= 10.0) * (0.8+pt*0.01) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0) * (0.85) + \
	  (abs(eta) > 4.0) * (0.00)

  }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  # tracking efficiency formula for muons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
	  (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 1.00) + \
	  (abs(eta) <= 1.2) * (pt > 1.0) * (1.00) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.8) * (pt > 0.2 && pt <= 1.0) * (pt*1.00) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.8) * (pt > 1.0) * (1.00) + \
	  (abs(eta) > 2.8 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + \
	  (abs(eta) > 2.8 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + \
	  (abs(eta) > 4.0) * (0.00)

  }
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  ## hadrons after having applied the tracking efficiency
  set InputArray  ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  # resolution formula for charged hadrons ,

  # from http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source trackMomentumResolution.tcl
}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons

  # taking something flat in energy for now, ECAL will take over at high energy anyway.
  # inferred from hep-ex/1306.2016 and 1502.02701
  set ResolutionFormula {

                        (abs(eta) <= 1.5)  * (energy*0.028) +
    (abs(eta) > 1.5  && abs(eta) <= 1.75)  * (energy*0.037) +
    (abs(eta) > 1.75  && abs(eta) <= 2.15) * (energy*0.038) +
    (abs(eta) > 2.15  && abs(eta) <= 3.00) * (energy*0.044) +
    (abs(eta) > 3.00  && abs(eta) <= 4.00) * (energy*0.10)}

}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  # resolution formula for muons

  # up to |eta| < 2.8 take measurement from tracking + muon chambers
  # for |eta| > 2.8 and pT < 5.0 take measurement from tracking alone taken from
  # http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source muonMomentumResolution.tcl
}


##############
# Track merger
##############

module Merger TrackMerger {
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}

########################################
#   Smear tracks
########################################

module TrackSmearing TrackSmearing {
  set InputArray TrackMerger/tracks
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray tracks
  set ApplyToPileUp true

  # from http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source trackResolutionCMS.tcl
}

########################################
#   Time Smearing
########################################

module TimeSmearing TimeSmearing {
  set InputArray TrackSmearing/tracks
  set OutputArray tracks

  # in sec
  set TimeResolution 30E-12
}

##################################
# Primary vertex reconstruction in 3D
##################################


module VertexFinder VertexFinder {
  set InputArray TimeSmearing/tracks

  set OutputArray tracks
  set VertexOutputArray vertices

  set MinPT 1.0
  set MinNDF 4
  set SeedMinPT 1.0

  set Sigma     3.0
  set MaxEta    10.0
  set GrowSeeds 1
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
  # add Branch InputArray BranchName BranchClass

  add Branch RHadronMerger/stableParticles stableGenParticle GenParticle
  add Branch TimeSmearing/tracks Track Track

  # add Branch VertexFinderDA4D/vertices Vertex4D Vertex
  add Branch VertexFinder/vertices Vertex Vertex
}
