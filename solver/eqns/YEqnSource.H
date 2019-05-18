// Add source terms for different precipitation models
if(precipitationModelName == "Oncul")
{
   volScalarField Sa(speciesThermo.Sa());
   YiEqn -= 
        fvm::Su
        (
            n * 3 * kv * Cg * rhod * (Wi/Ws) * m2 * pow((speciesThermo.Sa() - 1), Ng),
            Yi
        );
   
   Info<< "maxSource("<<Yi.name()<<") = " 
       << n * max(3 * kv * Cg * rhod * (Wi/Ws) * m2 * pow((speciesThermo.Sa() - 1), Ng)).value() << endl;

} else if(precipitationModelName == "Galbraith")
{
   //volScalarField SI(speciesThermo.SI());
   YiEqn -= 
       fvm::Su
       (
            n * 3 * ka * Cg * m2 * rhod * (Wi/Ws) * pow(speciesThermo.SI(), Ng),
            Yi
       );

   Info<< "maxSource("<<Yi.name()<<") = "
       << n * max(3 * ka * Cg * rhod *(Wi/Ws) * m2 * pow(speciesThermo.SI(), Ng)).value() << endl;

} else if(precipitationModelName == "Ye")
{
   //volScalarField S(speciesThermo.S());
   YiEqn -= 
       fvm::Su
       (
            n * 0.5 * ka * Cg * rhod * (Wi/Ws) * m2 * pow(speciesThermo.S(), Ng),
            Yi
       );

   Info<< "maxSource("<<Yi.name()<<") = "
       << n * max(0.5 * ka * Cg * rhod * (Wi/Ws) * m2 * pow(speciesThermo.S(), Ng)).value() << endl;

} else if(precipitationModelName == "Hanhoun")
{
   //volScalarField S(speciesThermo.S());
   YiEqn -= 
       fvm::Su
       (
            n * 3 *  kv * Cg * rhod * (Wi/Ws) * m2 * pow(speciesThermo.S(), Ng),
            Yi
       );

   Info<< "maxSource("<<Yi.name()<<") = "
       << n * max(3 *  kv * Cg * rhod * (Wi/Ws) * m2 * pow(speciesThermo.S(), Ng)).value() << endl;

} else 
{
    // Do nothing
}

