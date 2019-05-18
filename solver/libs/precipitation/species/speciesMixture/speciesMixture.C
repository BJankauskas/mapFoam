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

#include "speciesMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::speciesMixture::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < ROOTVSMALL)
    {
        FatalErrorIn
        (
            "void Foam::speciesMixture::"
            "correctMassFractions()"
        )
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesMixture::speciesMixture
(
    const dictionary& speciesThermoDict,
    const wordList& specieNames,
    const fvMesh& mesh
)
:
    species_(specieNames),
    Y_(species_.size()),
    W_(this->Y().size()),
    Z_(this->Y().size())
{
    Info<<"DEEP CODE IN SPECIES"<<endl;

    forAll(species_, i)
    {
        IOobject header
        (
            species_[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // Check if field exists and can be read
        if (header.headerOk())
        {
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField Ydefault
            (
                IOobject
                (
                    "Ydefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault
                )
            );
        }
    }

    correctMassFractions();

    forAll(species_, i)
    {
        const word& specieName = species_[i];
       // W_[i](speciesThermoDict.subDict(specieName).lookup("molWeight"));
        Z_[i] = speciesThermoDict.subDict(specieName).lookupOrDefault<scalar>("ionCharge", 0);

        W_.set
        (
            i,
            new dimensionedScalar
            (
                dimensionedScalar::lookupOrDefault
                (
                    "molWeight",
                    speciesThermoDict.subDict(specieName),
                    dimMass/dimMoles,
                    VSMALL
                )
            )
        );
    }
}

Foam::speciesMixture::speciesMixture
(
    const dictionary& speciesThermoDict,
    const fvMesh& mesh
)
:
    species_(speciesThermoDict.lookup("species")),
    Y_(species_.size()),
    W_(this->Y().size()),
    Z_(this->Y().size())
{
    forAll(species_, i)
    {
        IOobject header
        (
            species_[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // Check if field exists and can be read
        if (header.headerOk())
        {
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField Ydefault
            (
                IOobject
                (
                    "Ydefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault
                )
            );
        }
    }

    correctMassFractions();

    forAll(species_, i)
    {
        const word& specieName = species_[i];
        Z_[i] = speciesThermoDict.subDict(specieName).lookupOrDefault<scalar>("ionCharge", 0);

        W_.set
        (
            i,
            new dimensionedScalar 
            (
                dimensionedScalar::lookupOrDefault
                (
                    "molWeight",
                    speciesThermoDict.subDict(specieName),
                    dimMass/dimMoles,
                    VSMALL
                )
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesMixture::read
(
    const dictionary& speciesThermoDict
)
{}


// ************************************************************************* //
