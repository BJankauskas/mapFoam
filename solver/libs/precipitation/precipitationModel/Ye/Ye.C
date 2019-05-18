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

#include "Ye.H"
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

defineTypeNameAndDebug(Ye, 0);
addToRunTimeSelectionTable(precipitationModel, Ye, precipitationModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Ye::Ye
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
:
    precipitationModel(U, rho, phi),
    Cg_(this->populationBalanceProperties().subDict("univariateCoeffs")
                                           .subDict("growthModel")
                                           .lookup("Cg")),
    Ng_(readScalar(this->populationBalanceProperties().subDict("univariateCoeffs")
                                                      .subDict("growthModel")
                                                      .lookup("exponent"))),
    ka_(readScalar(this->speciesDict().lookup("ka"))),
    precipitant_(speciesDict().lookup("precipitant")),
    Ws_(this->speciesThermo().speciesComposition().W(precipitant_)),
    rhod_(mesh().lookupObject<volScalarField>("rhod"))
{
    Cg_.dimensions().reset(dimensionSet(0,1,-1,0,0,0,0));
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<Ye> Ye::New
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
{
    return autoPtr<Ye>
    (
        new Ye(U, rho, phi)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::Ye::precipitationSource
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

    volScalarField& Si = tSi.ref();

    // Species molecular weights
    const dimensionedScalar& Wi = 
        this->speciesThermo().speciesComposition().W(Y.name()); 

    // Second moment of CSD
    const volScalarField& m2 = mesh().lookupObject<volScalarField>("moment.2.populationBalance");

    Si = 0.5 * ka_ * Cg_ * rhod_ * (Wi/Ws_) * m2 * pow(this->speciesThermo().S(), Ng_);

    Info<< "maxSource("<<Y.name()<<") = " << max(Si).value() << endl;

    return tSi;
}

Foam::tmp<Foam::volScalarField>
Foam::Ye::alphaPrecipitationSource() const
{
    tmp<volScalarField> tAlphai 
    (
        new volScalarField
        (
            IOobject
            (
                "Alphai",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("zero", dimless/dimTime, 0)
        )
    );

    volScalarField& Alphai = tAlphai.ref();

    if(addAlphaSource_)
    {
        // Second moment of CSD
        const volScalarField& m2 = 
            mesh().lookupObject<volScalarField>("moment.2.populationBalance");

        Alphai = 0.5 * ka_ * Cg_ * m2 * pow(this->speciesThermo().S(), Ng_);

        Info<< "maxAlphaSource = " << max(Alphai).value() << endl;
    } else
    {
        // Do Nothing
    }

    return tAlphai;
}

bool Ye::read()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

