{
    precipitationModel->solvePBE();

    Switch calcSaueter = pbeDict.lookupOrDefault<Switch>("calculateSaueter", true);

    //TODO FUDGE
    const volScalarField& m0
    ( 
         mesh.lookupObject<volScalarField>("moment.0.populationBalance") 
    );
    const volScalarField& m1
    ( 
         mesh.lookupObject<volScalarField>("moment.1.populationBalance") 
    );
    const volScalarField& m2
    ( 
         mesh.lookupObject<volScalarField>("moment.2.populationBalance") 
    );
    const volScalarField& m3
    ( 
         mesh.lookupObject<volScalarField>("moment.3.populationBalance") 
    );
    const volScalarField& m4
    (
         mesh.lookupObject<volScalarField>("moment.4.populationBalance") 
    );

    alphaRatio == alpha1/(alpha1+ VSMALL);

//    Info<<"alphaRatio min = "<< min(alphaRatio).value() << "\n"
//        <<" alphaRatio max = "<< max(alphaRatio).value() << endl;
    
    forAll(alphaRatio, celli)
    {
        if(alpha1[celli] < 1e-6)
        {
            alphaRatio[celli] = 0.0;
        } 
    }

    if(calcSaueter)
    {
     
        const dimensionedScalar VS("VS", dimensionSet(0,0,0,0,0,0,0), 1e-9);
        const dimensionedScalar VS2("VS", dimensionSet(0,1,0,0,0,0,0), 1e-16);
        const dimensionedScalar VS3("VS", dimensionSet(0,0,0,0,0,0,0), 1e-16);
     
        //m4.dimensions().reset(dimless);
      //  m3.max(1e-9);
      //  m4.max(1e-9);

        forAll(m3, cellI)
        {
            if(m3[cellI] > SMALL)
            {
                d43[cellI] = m4[cellI] / m3[cellI];
            } else
            {
                d43[cellI] = 0;
            }
        }

        d32 = m4;
        d32 /= (m3 + VS3);
 
        Info<< "min d43 b4 = " << min(d32).value() 
            << " max b4 = " << max(d32).value() << endl;

        d32 = d32 * alphaRatio;

        Info<< "min d43 aft = " << min(d43).value() 
            << " max aft = " << max(d43).value() << endl;


      //  d32 == m3;
      //  d32 /= (m2 + VS2);
      //  
      //  d32 *= d32 * alphaRatio;
     
        Info<< "\nmin(m0) = "<< min(m0).value() << " max(m0) = " << max(m0).value() << "\n"
            << "min(m1) = "<< min(m1).value() << " max(m1) = " << max(m1).value() << "\n"
            << "min(m2) = "<< min(m2).value() << " max(m2) = " << max(m2).value() << "\n"
            << "min(m3) = "<< min(m3).value() << " max(m3) = " << max(m3).value() << "\n"
            << "min(m4) = "<< min(m4).value() << " max(m4) = " << max(m4).value() << "\n"
            << "min(VTTTT) = "<< min(g*d43*d43*(mixture.rhod() - mixture.rhoc())/(18 * mixture.mu())).value() 
            << " max(VTTTT) = " << max(g*d43*d43*(mixture.rhod() - mixture.rhoc())/(18 * mixture.mu())).value() << "\n"
            << "Volume averaged d43: "<< d43.weightedAverage(mesh.Vsc()).value()  << "\n" << endl;
            //<< "Volume averaged d32: "<< d32.weightedAverage(mesh.Vsc()).value()  << "\n"<< endl;
 
     
      //  scalar volAverage = 0;
      //  scalar totalVolume =0;
     
      //  forAll(d43, cellI)
      //  {
      //       totalVolume += mesh.V()[cellI];
      //  }
     
      //  forAll(d43, cellI)
      //  {
      //       volAverage += d43[cellI]*mesh.V()[cellI]/totalVolume;
      //  }
     
      //  Info<< "Volume averaged normalized d43: " << volAverage/(2.0e-6) << "\n"
     
      Vt == g*d43*d43*(mixture.rhod() - mixture.rhoc())/(18 * mixture.mu());

      Info<< "min(Vt) = "<< min(Vt).value() << " max(Vt) = "<< max(Vt).value() << endl;

    }

}
