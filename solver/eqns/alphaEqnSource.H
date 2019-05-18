// Add source terms for different precipitation models
if(alphaPrecipitationModel == "Oncul")
{
   alpha1Eqn -= 
        fvm::Su
        (
            3 * kv * Cg * m2 * pow((speciesThermo.Sa() - 1), Ng),
            alpha1 
        );
   
   Info<< "maxSource(alpha) = "
       << max(3 * kv * Cg * m2 * pow((speciesThermo.Sa() - 1), Ng)).value() << endl;

} else if(alphaPrecipitationModel == "Galbraith")
{
   alpha1Eqn -= 
       fvm::Su
       (
            3 * ka * Cg * m2 * pow(speciesThermo.SI(), Ng),
            alpha1
       );

   Info<< "maxSource(alpha) = "
       << max(3 * ka * Cg * m2 * pow(speciesThermo.SI(), Ng)).value() << endl;

} else if(alphaPrecipitationModel == "Ye")
{
   alpha1Eqn -= 
       fvm::Su
       (
            0.5 * ka * Cg * m2 * pow(speciesThermo.S(), Ng),
            alpha1
       );

   Info<< "maxSource(alpha) = "
       << max(0.5 * ka * Cg * m2 * pow(speciesThermo.S(), Ng)).value() << endl;

} else if(alphaPrecipitationModel == "Hanhoun")
{
   alpha1Eqn -= 
       fvm::Su
       (
            3 *  kv * Cg * m2 * pow(speciesThermo.S(), Ng),
            alpha1
       );

   Info<< "maxSource(alpha) =  "
       << max(3 *  kv * Cg * m2 * pow(speciesThermo.S(), Ng)).value() << endl;

} else 
{
    // Do nothing
}

