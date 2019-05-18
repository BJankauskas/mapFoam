/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "tFudge.H"
#include "addToRunTimeSelectionTable.H"

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(tFudge, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, tFudge, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::tFudge::tFudge
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    a_("a", dimless, dict),
    a1_("a1", dimless, dict),
    V0_("V0", dimVelocity, dict),
    residualAlpha_("residualAlpha", dimless, dict),
    Vt_(this->mixture().U()),
    g_("g", dimAcceleration, dict) //TODO FUDGE - should just used regular g field; sort it out 
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::tFudge::~tFudge()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::tFudge::correct()
{

    const volScalarField& dd = this->mixture().d32();
    //const volScalarField& dd = this->mixture().d43();
    const volScalarField& mu = this->mixture().mu();

    Udm_ =
        (rhoc_/rho())
       *(g_*dd*dd*(rhod_-rhoc_)/(18*mu))
       *(
            exp(-a_*max(alphad_ - residualAlpha_, scalar(0)))
          - exp(-a1_*max(alphad_ - residualAlpha_, scalar(0)))
        );

    Info<< "min(Udm) = " << min(Udm_).value() << " max(Udm) = " << max(Udm_).value() <<"\n"
        << "max(Vt) = "<< max((g_*dd*dd*(rhod_ - rhoc_)/(18 * mu))).value() << "\n"
        << "max(rhox/rho()) = "<< max(rhoc_/rho()).value() << "\n"
        << "min(mu) = "<< min(mu).value() << " max(mu) = "<< max(mu).value() << endl;

}


// ************************************************************************* //
