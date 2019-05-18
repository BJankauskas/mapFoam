tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        rhoPhi,
        mesh.divScheme("div(phi,Yi)")
    )
);

{
//    precipitation->correct(); 
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    
    // Update saturation indices
    SI == speciesThermo.SI();
    S == speciesThermo.S();
    Sa == speciesThermo.Sa();

    scalar n;

    forAll(Y, i)
    {

        const dimensionedScalar& Wi = composition.W(Y[i].name());
       
        if (Y[i].name() != inertSpecie)
        {
            
            volScalarField& Yi = Y[i];

            // Source term sign depending on species being solved
            if(Yi.name() == precipitationModel->precipitant())
            { n = 1; } else { n = -1; }

            // Declare RHS of YEqn
            fvScalarMatrix YiEqn
            (
                  fvm::ddt(rho, Yi)
                + mvConvection->fvmDiv(rhoPhi, Yi)
                - fvm::laplacian(turbulence->muEff(), Yi)
                ==
                  n * fvm::Su(precipitationModel->precipitationSource(Yi), Yi)
            );

      //      #include "YEqnSource.H"

            YiEqn.relax();

            YiEqn.solve(mesh.solver("Yi"));

            Yi.max(0.0);
            Yt += Yi;

        }
        else
        {
            inertIndex = i;
        }
    }
     
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

    speciesThermo.printThermoParams();
    speciesThermo.updateConcentrations();
}
