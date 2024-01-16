# Planetary fluvial sediment transport model 
Contact: Lisanne Braat (lisannebraat@gmail.com)
Braat et al. 2024 - Gravity-Driven Differences in Fluvial Sediment Transport on Mars and Earth

This MATLAB model accompanies the paper "Gravity-Driven Differences in Fluvial Sediment Transport Fluxes on Mars and Earth" by Lisanne Braat, Muriel Z. M. Brückner, Elliot Sefton-Nash and Michael P. Lamb. The model was written in MATLAB version R2022b.
The model was used to compare fluvial sediment transport on Mars and Earth and to specifically analyse the effect of gravity.
The boundary/input conditions (independent variables) can be easily modified to test other conditions and several sediment transport relations are included for comparison. 
Please refer to the paper of more information.

The model consists of the following scripts:
- FlowOnMars1.m
- FlowOnMars2_InitiationMotion.m
- FlowOnMars3_Sediment.m
- FlowOnMars4_Transport_h.m
- FlowOnMars4_Transport_Q.m
- FlowOnMars5_GrainSizeDistribution.m
- FlowOnMars6_TransportMixture.m
- FlowOnMars7_IndependentVariables.m

The name of the script related to the output highest level output that is generated in the form of figures. The scripts can be run independently, except that FlowOnMars5 is required to run before running FlowOnMars6 or FlowOnMars7. The scripts are described in more detail below.

-------------------------------------

FlowOnMars1.m
Based on a set of input parameters (channel width, slope, water density, gravity, temperature, roughness height and water depth or discharge), some hydrological parameters are calculated and visualised (water depth, hydraulic radius, roughness, velocity, bed shear stress, shear velocity, Froude number, Reynolds number, laminar sublayer thickness).

FlowOnMars2_InitiationMotion.m
This script compares a large set of relations for the initiation of sediment motion.

FlowOnMars3_Sediment.m
Based on the same input as for FlowOnMars1 + sediment density and a grain size range some parameters are calculated in relation to sediment (settling velocity, particle Reynolds number, Bonnefile number, movability number, advection length, Shields number, critical Shields number, critical shear stress, critical velocity, critical shear velocity)

FlowOnMars4_Transport_h.m
Input like in FlowOnMars1 + FlowOnMars3 is used to calculate bedload transport, suspended load transport and total transport fora specified grain size range. Many sediment transport equations are implemented to choose from or to compare. *_h required a water depth input.

FlowOnMars4_Transport_Q.m
Input like in FlowOnMars1 + FlowOnMars3 is used to calculate bedload transport, suspended load transport and total transport for a specified grain size range. Many sediment transport equations are implemented to choose from or to compare. *_Q required a water discharge input.

FlowOnMars5_GrainSizeDistribution.m
Generates a lognormal sediment distribution for FlowOnMars6 and FlowOnMars7. You can generate your own distribution as long as it is in the same format as this output.

FlowOnMars6_TransportMixture.m
Calculated sediment transport based on a sediment mixture. This script only includes one method of calculating sediment transport.

FlowOnMars7_IndependentVariables.m
Calculated sediment transport based on a sediment mixture but for various independent parameters (water depth, discharge, velocity, and shear stress), so multiple hydrodynamic input parameters can be compared. This script only includes one method of calculating sediment transport.

