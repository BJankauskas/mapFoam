/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alberto Passalacqua
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

#include "powerLawNucleation.H"
#include "addToRunTimeSelectionTable.H"
//#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(powerLawNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        powerLawNucleation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::powerLawNucleation::powerLawNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    coeffsDict_(dict.subDict(this->type() + "Coeffs")),
    Kg_(coeffsDict_.lookup("nucleationRate")),
    exponent_(readScalar(coeffsDict_.lookup("exponent"))),
    J_
    (
        IOobject
        (
            "nucleationSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/(dimVolume * dimTime), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::powerLawNucleation::~powerLawNucleation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::populationBalanceSubModels::nucleationModels::powerLawNucleation
::tInduction(const volUnivariateMoment& moment) const
{
    Info<<"Induction Time: min(T_induction) = " << 1/(max(J_).value() + VSMALL)
        <<" max(T_induction) = " << 1/(min(J_).value() + VSMALL) << "\n" << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::powerLawNucleation
::nucleationSource(const volUnivariateMoment& moment) //const
{

    const volScalarField& SI_ = mesh_.lookupObject<volScalarField>("SI");

    tmp<volScalarField> tJ
    (
        new volScalarField
        (
            IOobject
            (
               "J", 
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", moment.dimensions()/dimTime, 1.0)
        )
    );

    volScalarField& J = tJ.ref();

    scalar abscissaNucleation = 0;

    J *= pow(SI_, exponent_);
    if(moment.order() == 0)
    {   
        J_ == Kg_ * J;

        Info << "min(J) = " << min(Kg_*J).value()
             << " max(J) = " << max(Kg_*J).value() << endl;

        this->tInduction(moment);
    }
    
    return Kg_ * pow(abscissaNucleation, moment.order()) * tJ;
}

// ************************************************************************* //

