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

#include "OnculGrowth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(OnculGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        OnculGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::OnculGrowth
::OnculGrowth
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

Foam::populationBalanceSubModels::growthModels::OnculGrowth
::~OnculGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::OnculGrowth::Kg
(
    const volScalarField& abscissa
) const
{

    const fvMesh& mesh_ = abscissa.mesh();
    const volScalarField& Sa_ = mesh_.lookupObject<volScalarField>("Sa"); 


    // TODO FUDGE -> this is a nasty way of doing it
    volScalarField G(Sa_);
    G = pow((G - 1), exponent_);

    dimensionedScalar oneAbs
    (
        "oneAbs",
        dimVolume/sqr(abscissa.dimensions()),
        1.0
    );

    return Cg_*pos(-abscissa + maxAbscissa_)
        *pos(abscissa - minAbscissa_)*oneAbs*G;
}

// ************************************************************************* //
