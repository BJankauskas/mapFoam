/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "precipitate.H"
#include "Time.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(precipitate, 0);
addToRunTimeSelectionTable(precipitationModel, precipitate, precipitationModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

precipitate::precipitate
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
:
    precipitationModel(U, rho, phi)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<precipitate> precipitate::New
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
{
    return autoPtr<precipitate>
    (
        new precipitate(U, rho, phi)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::precipitate::precipitationSource
(
    volScalarField& Y
) const
{
    tmp<volScalarField> tSi
    (
        new volScalarField
        (
            IOobject
            (
                "Si",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

//    volScalarField& Si = tSi.ref();
//
//   //  word precipitationModelName_ = this->lookup("precipitationModelName");
//
//    // Growth rate and exponent
//    dimensionedScalar Cg = this->populationBalanceProperties().subDict("univariateCoeffs")
//                                                               .subDict("growthModel")
//                                                               .lookup("Cg");
//    const scalar& Ng = readScalar(this->populationBalanceProperties().subDict("univariateCoeffs")
//                                                               .subDict("growthModel")
//                                                               .lookup("exponent"));
//  
//    Cg.dimensions().reset(dimensionSet(0,1,-1,0,0,0,0));
//
//
//    Info<< "b4 ks" << en>populationBalanceProperties().subDict("univariateCoeffs")
//                                                               .subDict("growthModel")
//                                                               .lookup("Cg");
//    const scalar& Ng = readScalar(this->populationBalanceProperties().subDict("univariateCoeffs")
//                                                               .subDict("growthModel")
//                                                               .lookup("exponent"));
//  
//    Cg.dimensions().reset(dimensionSet(0,1,-1,0,0,0,0));
//
//
//    Info<< "b4 ks" << endl;
//
//    // Shape factors
//    scalar kv = readScalar(speciesDict().lookup("kv"));
//    scalar ka = readScalar(speciesDict().lookup("ka"));
//
//
//    // Shape factors
//    scalar kv = readScalar(speciesDict().lookup("kv"));
//    scalar ka = readScalar(speciesDict().lookup("ka"));
//
//    Info<< "b4 precipitatnt" << endl;
//    // Precipitant name
//    const word& precipitant = speciesDict().lookup("precipitant");
//
//    // Species molecular weights
//    const dimensionedScalar& Wi = this->speciesThermo().speciesComposition().W(Y.name()); 
//    const dimensionedScalar& Ws = this->speciesThermo().speciesComposition().W(precipitant);
//
//    // Crystal density
//    const volScalarField& rhod = mesh().lookupObject<volScalarField>("rhod");
//    //const dimensionedScalar& rhod(Wi); 
//        //= mesh().lookupObject<dimensionedScalar>("rhod");
//
//    Info<< "b4 m2" << endl;
//    // Second moment of CSD
//    const volScalarField& m2 = mesh().lookupObject<volScalarField>("moment.2.populationBalance");
//
//    // Temporary solution TODO
//    volScalarField Sa(this->speciesThermo().Sa());
//    volScalarField S(this->speciesThermo().S());
//    volScalarField SI(this->speciesThermo().SI());
//
//
//    Info<< "after Sa" << endl;
//    // +/- factor
//    scalar n;
//    if(Y.name() == precipitant)
//    { n = 1; } else { n = -1; }

  //  if(precipitationModelName_ == "Oncul")
  //  {
  //     //speciesSource_ == n * 3 * kv * Cg * rhod * (Wi/Ws) * m2 * pow((speciesThermo().Sa() - 1), Ng);
  //     Si == n * 3 * kv * Cg * rhod * (Wi/Ws) * m2 * pow((Sa - 1), Ng);
  //  } else if(precipitationModelName_ == "Galbraith")
  //  {
  //    // speciesSource_ == n *0.5 * ka * Cg * rho() * pow(speciesThermo().SI(), Ng);
  //     Si == n *0.5 * ka * Cg * this->speciesThermo().rho() * pow(SI, Ng);
  //  } else if(precipitationModelName_ == "Ye")
  //  {
  //     //speciesSource_ == n * 0.5 * ka * Cg * m2 * pow(speciesThermo().S(), Ng);
  //     Si == n * 0.5 * ka * Cg * m2 * pow(S, Ng);
  //  } else if(precipitationModelName_ == "Hanhoun")
  //  {
  //     //speciesSource_ == n * 3 *  kv * Cg * rhod * (Wi/Ws) * m2 * pow((speciesThermo() - 1), Ng) 
  //     Si == n * 3 *  kv * Cg * rhod * (Wi/Ws) * m2 * pow((Sa - 1), Ng);
  //  } else 
  //  {
  //      // Do nothing
  //  }



    return tSi;
}

/*
Foam::tmp<Foam::fvScalarMatrix>
Foam::precipitate::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimMass/dimTime)
    );

    return tSu;
}
*/

//void Foam::precipitate::solvePBE()
//{
//    this->populationBalanceModel_->solve();
//}


bool precipitate::read()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

