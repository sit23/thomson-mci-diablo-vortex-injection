The original code was developed as a joint project in MAE 223 (CFD), taught by
T. Bewley, at USan Diego (spring 2001), and was provided to Stephen Thomson by 
Andy Thompson and Emma Boland in October 2011.

The present version of the code encorperates a vortex injection scheme developed by
Stephen Thomson and Michael McIntyre. Further information can be found on their websites:
http://www.damtp.cam.ac.uk/user/sit23/ 
http://www.atm.damtp.cam.ac.uk/mcintyre/
Any enquires should be sent to stephen.i.thomson@gmail.com
Formal documentation discussing the vortex injection scheme will be available in due course,
although details are discussed in our submitted paper, available on our websites.

****************************************************************************************
This code is free software; you can redistribute it and/or modify it
under the terms of the GNU General PubliLicense as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version. This code is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General PubliLicense for more details. You should have received a
copy of the GNU General PubliLicense along with this code; if not,
write to the Free Software Foundation, Inc., 59 Temple Place - Suite
330, Boston, MA 02111-1307, USA.
****************************************************************************************

More comprehensive documentation will be provided in due course, but the code is displayed in its current state for the interested reader. Many IDL routines have been written by Stephen Thomson to analyse the code's output, and to convert model units to physical units. Please email Stephen on stephen.i.thomson@gmail.com if you would like to see them.

****************************************************************************************
UNIT CONVERSION

The conversions from model units to physical units is as follows:

10,000km is equal to 32 /pi model units. 
0.0005 is 50 seconds. So 1 second is 10^-5 model time units.

This conversion means that model speeds are approximately equal to values in ms^-1.
****************************************************************************************

INPUT VARIABLES

	BETA = 0.00 Beta value for particular latitude and planet. 
	KAPPA = 0.00 Leave unchanged to reproduce results
	NU = 0.00 Leave unchanged to reproduce results
    HX = 0.00 Constant slope in x direction for topography.
     HY = 0.0 Constant slope in y direction for topography.
	LX = 201.062 X domain length
	LY = 100.531 Y domain length
	CSX = 1.75 Leave unchanged to reproduce results
	CSY = 1.75 Leave unchanged to reproduce results
	N_TIME_STEPS = 40000000
	DELTA_T = 0.0005
	RESET_TIME = .FALSE.
	NUM_PER_DIR = 2 Leave unchanged
	TIME_AD_METH = 1 Leave unchanged
	VERBOSITY = 1 Leave unchanged
	SAVE_FLOW_INT = 50000 How often to output FFT(psi) for the purposes of restarting the code. (in timesteps)
	SAVE_PHYS_INT = 50000 How often to output the physical fields (PV, psi, etc)
      SAVE_PROF_INT = 50000 How often to output profiles (zonal means etc)
	SAVE_STATS_INT = 50000 How often to output energies
	CREATE_NEW_FLOW = .TRUE.
	IC_TYPE = 1 Leave unchanged
	ENERGY = 0.5 Makes no difference with current initial condition set up
	IC_1 = 0.0 Leave unchanged
	IC_2 = 0.0 Leave unchanged
	BC1_TYPE = 1 Leave unchanged
	BC2_TYPE = 1 Leave unchanged
    HAMP = 2.1875 !Part of topography amplitude.
    KCUT = 85 !Cutoff wavenumber in Smith filter. 
    FILTER_EXP = 8 !filter exponent in Smith filter
	LDSQD = 145.53 !ld^2 in model units
	KT = 0.0625 !wavevector of psi_2 variation in y direction
	SR = 2.0 !constant determining injection strength. 
	CMC = 0.000007138217212 !Second constant determining injection strength. Always appears in ratio sr/cmc. 
	a = 1.0 Leave unchanged
	b = 1.0 Leave unchanged
	to = 0.001 Time overwhich injections take place (make it multiple of 2*delta_t for computational-mode avoidance.)
	xo = 1.570 Injected storm radius
	psilim = -11569.92 Psilim value (note model definition slightly differnt to paper)
	tcool = 12.6752  Irrelevant
	bias_factor = 16.0 1/b_max
	time_bet_storms = 1728.0 !Max number of timesteps between injections -4. (real time gap = 4.+time_bet_storms*rand())
	q_thresh=35.0 (max q' strength of vortex in model units).
	use_new_storms = .TRUE. / !Make this false if you want to read in old storm list. 
****************************************************************************************
A FEW OF THE OUTPUT VARIABLES

fort.12 psi
fort.14 PV without beta*y+psi_2/ldsqd
fort.16 PV with beta*y+psi_2/ldsqd
fort.90 A list of all of the injected storms (2 records per storm, as each injection takes two timesteps)
****************************************************************************************

TIME ADVANCEMENT

The equation the model solves is:

D(q)/Dt = F(x,y,t)

where q=\nabla^2 \psi + \beta y -(\psi-\psi_deep)/ldsqd. 
and F(x,y,t) is the forcing due to vortex injections.

The time advancement advances the relative vorticity in spectral space, which (with the appropriate scaling applied in the model) is exactly equivalent to timestepping the full PV.

To do the time advancement, a field CR is calculated, which is the RHS of the equation \partial(FT(q_relative))/\partial t = CR

CR is therefore a function of x and y wavenumbers (k_x,k_y). So, CR(k_x,k_y)=FT(F - u.grad q)/(1+(1/(ldsqd*(k_x^2 + k_y^2)))). (F- u.grad q) is calculated in physical space, and subsequently fourier transformed. The factor 1/(1+(1/(ldsqd*(k_x^2 + k_y^2)))), if taken to the LHS would turn FT(q_rel) into FT(q_rel - psi/ldsqd). As there is a time derivative applied to this, one could make it FT(q_rel + beta y -(psi-psi_deep)/ldsqd), which would be the required timestepping of the full PV.

The time advancement is a simple leap-frog, with weak Robert filter. The filter is only applied when there is NOT an injection currently taking place, to avoid exciting a computational mode.
****************************************************************************************

Stephen Thomson (December 2014)
stephen.i.thomson@gmail.com



