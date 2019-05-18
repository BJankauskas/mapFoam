# mapFoam 

This solver is tailored to modelling of precipitation reactions of Struvite (magnesium ammonium phosphate).
The solver is aimed to capture the hydrodynamics, equilibrium reaction thermodynamics and the evolution of the newly formed dispersed phase.

In order to run the solver there are a few prerequisites.

1) Written under OpenFOAM-4.x version. Later versions would require some porting although it shouldn't be too much, other than some minor syntaxial changes
2) OpenQBMM-2.0 framework was used for Population balance equation solution. The group actively develops the framework, however the later versions have had some major
   changes done, hence the version of the framework has been left as it was at the start of the solver development.
   - Please use the OpenQBMM package in this repository as there has been some minor changes done to couple of classes; nothing major, but few classes required change
     in access permisions, so that the solver would be able to access certain terms when required. This is a terrible design, but it was left as it is for the time being.



Notes: 
    - Solver is provided mainly as proof of concept, given that no data was available to perform validation studies (see my thesis). 
    - The aim of the solver is to build a framework for precipitation processes, with the focus on struvite, however the principles would remain the same for different reactive processes.
      If the solver is used for different chemistry, then thermodynamics classes have to be rewritten for that particular system. Unfortunately, I don't think there is a way to write it in
      a generic way, because equilibrium speciation of each system would be different and therefore cannot be generalised;
    - Solver is "designed" in an object oriented manner as much as possible to make it flexible for future adjustments, nonetheless most likely there are certain decisions made that would
      make the user roll their eyes, therefore again I would like to reiterate that the aim of the solver is to be used as a skeleton for the framework which can and should be dismantled
      and then pulled back into something that is useful for a new purpose.

TODO:
    - Fix the initial condition calculations. IDEA is to specify required concentrations for the chemical species in one of the dictionaries, which then used at the required inlets. (NOT SURE HOW YET);
    - d43 clipping; the blow up fix is currently a mystery; current hack gives jagged results
    - Turn the OpenQBMM realisability warnings and figure out how reaction algorithm can and should be adjusted so that the moments inversion is more stable.
    - Bingham plastic is probably a wrong choice for dilute phase; Is it not just a Newtonian? (This is the limitation of driftFluxFoam as it was written for sludge and here we don't really have that;
                                                                                                check with literature)

