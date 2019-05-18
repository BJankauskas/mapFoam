/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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

#include "GalbraithAggregation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(GalbraithAggregation, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        GalbraithAggregation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::GalbraithAggregation
::GalbraithAggregation
(
    const dictionary& dict
)
:
    aggregationKernel(dict),
    exponent_(readScalar(dict.lookup("exponent")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::GalbraithAggregation
::~GalbraithAggregation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::aggregationKernels::GalbraithAggregation::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{   
    const fvMesh& mesh_ = abscissa1.mesh();
    const volScalarField& SI_ = mesh_.lookupObject<volScalarField>("SI");   // Basing on Galbraith et al.

    tmp<volScalarField> tAgg
    (
        new volScalarField
        (
            IOobject
            (
                "GalbraithAggregationK",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            SI_
         )
    );

    volScalarField& Agg = tAgg.ref();

    Agg = pow(Agg, exponent_);
    Agg.dimensions().reset(pow3(abscissa1.dimensions())/dimTime);
       
   // Info<<"Agg = "<< min(Ca_.value() * Agg).value() << endl;

    return Ca_.value()*tAgg;
}

// ************************************************************************* //
