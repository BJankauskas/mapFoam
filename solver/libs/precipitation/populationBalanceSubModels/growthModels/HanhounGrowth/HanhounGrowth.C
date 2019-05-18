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

#include "HanhounGrowth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(HanhounGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        HanhounGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::HanhounGrowth
::HanhounGrowth
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

Foam::populationBalanceSubModels::growthModels::HanhounGrowth
::~HanhounGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::HanhounGrowth::Kg
(
    const volScalarField& abscissa
) const
{

    const fvMesh& mesh_ = abscissa.mesh();
    const volScalarField& S_ = mesh_.lookupObject<volScalarField>("S");     // Based on Hanhoun et al.

    volScalarField G(S_);
    G.dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
    G = pow(G, exponent_);

    dimensionedScalar oneAbs
    (
        "oneAbs",
        dimVolume/sqr(abscissa.dimensions()),
        1.0
    );

  //  Info<<"Cg = "<< Cg_.value() << "\n"
  //      <<"S = " << max(S_).value() << "\n"
  //      <<"G = "<< max(Cg_*G).value() << endl;

    return Cg_*pos(-abscissa + maxAbscissa_)
        *pos(abscissa - minAbscissa_)*oneAbs*G;
}

// ************************************************************************* //
