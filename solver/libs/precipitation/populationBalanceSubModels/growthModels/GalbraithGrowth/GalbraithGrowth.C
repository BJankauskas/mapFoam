/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Matteo Icardi
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "GalbraithGrowth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(GalbraithGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        GalbraithGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::GalbraithGrowth
::GalbraithGrowth
(
    const dictionary& dict
)
:
    growthModel(dict),
    minAbscissa_(dict.lookup("minAbscissa")),
    maxAbscissa_(dict.lookup("maxAbscissa")),
    exponent_(readScalar(dict.lookup("exponent")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::GalbraithGrowth
::~GalbraithGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::GalbraithGrowth::Kg
(
    const volScalarField& abscissa
) const
{

    const fvMesh& mesh_ = abscissa.mesh();
    const volScalarField& SI_ = mesh_.lookupObject<volScalarField>("SI");   // Basing on Galbraith et al.

    tmp<volScalarField> tG
    (
        new volScalarField
        (
            IOobject
            (
                "G",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            SI_
        )
    );

    volScalarField& G = tG.ref();

    G.dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
    G = pow(G, exponent_);

    dimensionedScalar oneAbs
    (
        "oneAbs",
        dimVolume/sqr(abscissa.dimensions()),
        1.0
    );

    return Cg_*pos(-abscissa + maxAbscissa_)
        *pos(abscissa - minAbscissa_)*oneAbs*tG;
}

// ************************************************************************* //
