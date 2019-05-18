/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "volFields.H"
#include "precipitationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(precipitationModel, 0);
    defineRunTimeSelectionTable(precipitationModel, precipitationModel);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::precipitationModel::precipitationModel
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
:
    IOdictionary 
    (
        IOobject
        (
            "precipitationProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(U.mesh()),
   // precipitationModelName_(this->lookup("precipitationModelName")),
   // alphaPrecipitationModelName_(this->lookupOrDefault<word>("alphaPrecipitation", "none")),
    speciesDict_(this->subDict("speciesMixture")),
    speciesThermoModel_(speciesThermoModel::New(U, rho, speciesDict_)),
    populationBalanceProperties_
    (
        IOobject
        (
            "populationBalanceProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    populationBalanceModel_
    (
        populationBalanceModel::New
        (
            "populationBalance", populationBalanceProperties_, U, phi
        )
    ),
    addAlphaSource_(this->lookupOrDefault<Switch>("addAlphaSource", true)),
    precipitant_(speciesDict_.lookup("precipitant"))
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<precipitationModel> precipitationModel::New
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "precipitationProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("simulationType")
    );

    Info<< "Selecting precipitation model type " << modelType << endl;

    precipitationModelConstructorTable::iterator cstrIter =
        precipitationModelConstructorTablePtr_->find(modelType);

    if (cstrIter == precipitationModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "precipitationModel::New(const volVectorField&, "
            "const surfaceScalarField&, precipitationModel&, const word&)"
        )   << "Unknown precipitationModel type "
            << modelType << nl << nl
            << "Valid precipitationModel types:" << endl
            << precipitationModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<precipitationModel>
    (
        cstrIter()(U, rho, phi)
    );
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::precipitationModel::~precipitationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::precipitationModel::solvePBE()
{
    this->populationBalanceModel_->solve();
}

bool Foam::precipitationModel::read()
{
    return true;
}


// ************************************************************************* //

} // End namespace Foam 

// ************************************************************************* //

