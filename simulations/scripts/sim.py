from crpropa import *
import numpy as np
import sys

RESULTS_DIR = '../results'

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZS = [1, 2, 7, 14, 56]

A = int(sys.argv[1])
Z = int(sys.argv[2])

for iZ in range(len(ZS)):
	if Z == ZS[iZ]:
		break

src_data = np.genfromtxt(f'starburst_galaxies.dat', dtype = None, encoding = None)

for isrc in range(32, len(src_data)):

	sim = ModuleList()
	sim.add( SimplePropagation(1*kpc, 10*Mpc) )
	sim.add( Redshift() )
	sim.add( PhotoPionProduction(CMB()) )
	sim.add( PhotoPionProduction(IRB_Saldana21()) )
	sim.add( PhotoDisintegration(CMB()) )
	sim.add( PhotoDisintegration(IRB_Saldana21()) )
	sim.add( NuclearDecay() )
	sim.add( ElectronPairProduction(CMB()) )
	sim.add( ElectronPairProduction(IRB_Saldana21()) )
	sim.add( MinimumEnergy(1 * EeV) )

	obs = Observer()
	obs.add( Observer1D() )
	if src_data[isrc][2] == 'NGC4038/9':
		output = TextOutput(f'{RESULTS_DIR}/{PARTICLES[iZ]}/events_NGC4038_9.txt', Output.Event1D)
	else:
		output = TextOutput(f'{RESULTS_DIR}/{PARTICLES[iZ]}/events_{src_data[isrc][2]}.txt', Output.Event1D)
	obs.onDetection( output )
	sim.add( obs )

	source = Source()
	source.add( SourcePosition(src_data[isrc][8] * Mpc) )
	source.add( SourceRedshift1D() )
	source.add( SourcePowerLawSpectrum(1 * EeV, 10000 * EeV, -1) )
	source.add( SourceParticleType(nucleusId(A, Z)) )

	sim.setShowProgress(True)
	sim.run(source, 100000, True)
	output.close()
