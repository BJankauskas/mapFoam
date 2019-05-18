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

#include "None.H"
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

defineTypeNameAndDebug(None, 0);
addToRunTimeSelectionTable(precipitationModel, None, precipitationModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

None::None
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
:
    precipitationModel(U, rho, phi)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<None> None::New
(
    volVectorField& U,
    volScalarField& rho,
    surfaceScalarField& phi
)
{
    return autoPtr<None>
    (
        new None(U, rho, phi)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::None::precipitationSource(volScalarField& Y) const
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

    return tSi;
}

Foam::tmp<Foam::volScalarField>
Foam::None::alphaPrecipitationSource() const
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

    return tAlphai;
}

bool None::read()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
