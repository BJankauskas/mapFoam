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
#include "speciesThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(speciesThermoModel, 0);
    defineRunTimeSelectionTable(speciesThermoModel, speciesThermoModel);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesThermoModel::speciesThermoModel
(
    volVectorField& U,
    volScalarField& rho,
    dictionary& speciesDict
)
:
    mesh_(U.mesh()),
    rho_(rho),
    speciesDict_(speciesDict),
    speciesMixture_(*new speciesMixture(speciesDict_, mesh_)), 
    printParams_(speciesDict_.lookupOrDefault<Switch>("printThermoParams", true))
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<speciesThermoModel> speciesThermoModel::New
(
    volVectorField& U,
    volScalarField& rho,
    dictionary& speciesDict
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
        ).lookup("speciesThermoModel")
    );

    Info<< "Selecting species transport model type " << modelType << endl;

    speciesThermoModelConstructorTable::iterator cstrIter =
        speciesThermoModelConstructorTablePtr_->find(modelType);

    if (cstrIter == speciesThermoModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "speciesThermoModel::New(const volVectorField&, "
            "const surfaceScalarField&, speciesThermoModel&, const word&)"
        )   << "Unknown speciesThermoModel type "
            << modelType << nl << nl
            << "Valid speciesThermoModel types:" << endl
            << speciesThermoModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<speciesThermoModel>
    (
        cstrIter()(U, rho, speciesDict)
    );
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::speciesThermoModel::~speciesThermoModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::speciesThermoModel::read()
{
    return true;
}


// ************************************************************************* //

} // End namespace Foam 

// ************************************************************************* //
