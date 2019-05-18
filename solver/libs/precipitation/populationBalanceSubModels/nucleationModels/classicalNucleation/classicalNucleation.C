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

#include "classicalNucleation.H"
#include "addToRunTimeSelectionTable.H"
//#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(classicalNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        classicalNucleation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::classicalNucleation::classicalNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    coeffsDict_(dict.subDict(this->type() + "Coeffs")),
    B1_(coeffsDict_.lookup("B1")),
    B2_(coeffsDict_.lookupOrDefault<dimensionedScalar>("B2", B1_)),
    A1_(readScalar(coeffsDict_.lookup("A1"))),
    A2_(coeffsDict_.lookupOrDefault<scalar>("A2", A1_)),
    SaiMin_(coeffsDict_.lookupOrDefault<scalar>("SaiMin", 1)),
    SaiMax_(coeffsDict_.lookupOrDefault<scalar>("SaiMax", 1000)),
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
{
    Info<< "SaiMin = " << SaiMin_
        << " SaiMax = " << SaiMax_  << "\n"
        << "A1 = " << A1_ 
        << " A2 = " << A2_ << "\n"
        << "B1 = " << B1_ 
        << " B2 = " << B2_ << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::classicalNucleation::~classicalNucleation()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::populationBalanceSubModels::nucleationModels::classicalNucleation
::tInduction(const volUnivariateMoment& moment) const
{
    Info<<"Induction Time: min(T_induction) = " << 1/(max(J_).value() + VSMALL)
        <<" max(T_induction) = " << 1/(min(J_).value() + VSMALL) << "\n" << endl;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::classicalNucleation
::nucleationSource(const volUnivariateMoment& moment) //const
{

    const volScalarField& Sa_ = mesh_.lookupObject<volScalarField>("Sa");

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
            dimensionedScalar("zero", moment.dimensions()/dimTime, 0)
        )
    );

    volScalarField& J = tJ.ref();

    scalar abscissaNucleation = 0;  

    //TODO There is a bug at the walls, where J is set to 0
    forAll(Sa_, cellI)
    {
        const scalar Sai_ = Sa_[cellI];

        if(Sai_ > SaiMin_ && Sai_ < SaiMax_) {
            J[cellI] = B1_.value() * exp(-A1_/(pow(log(Sai_), 2))); 
        } else if(Sai_ >= SaiMax_){
            J[cellI] = B2_.value() * exp(-A2_/(pow(log(Sai_), 2))); 
        } else {
            J[cellI] = 0;
        }
    }

    if(moment.order() == 0)
    {
        J_ == J;
        
        Info << "\n" << "Nucleation source: min(J) = " << min(J).value()
            << " max(J) = " << max(J).value() << endl;

        this->tInduction(moment);
    }

   
    return pow(abscissaNucleation, moment.order()) * tJ;
}

// ************************************************************************* //

