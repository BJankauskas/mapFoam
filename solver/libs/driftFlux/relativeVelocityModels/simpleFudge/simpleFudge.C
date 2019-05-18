/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "simpleFudge.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(simpleFudge, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, simpleFudge, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simpleFudge::simpleFudge
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    a_("a", dimless, dict),
    V0_("V0", dimVelocity, dict),
    residualAlpha_("residualAlpha", dimless, dict),
    Vt_(this->mixture().U()),      // TODO FUDGE
    g_("g", dimAcceleration, dict) // TODO FUDGE
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simpleFudge::~simpleFudge()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::simpleFudge::correct()
{

    //const volScalarField& dd = this->mixture().d32();
    const volScalarField& dd = this->mixture().d43();
    const volScalarField& mu = this->mixture().mu();

    Udm_ = 
        (rhoc_/rho())
        *(g_*dd*dd*(rhod_ - rhoc_)/(18 * mu))
        *pow(scalar(10), -a_*max(alphad_, scalar(0)));


  //  Info<< "min(Udm) = " << min(Udm_).value() << " max(Udm) = " << max(Udm_).value() <<"\n"
  //      << "mag(Vt) = "<< mag(min((g_*dd*dd*(rhod_ - rhoc_)/(18 * mu))).value()) << "\n"
  //      << "max(rhod - rhoc) = "<<(rhod_ - rhoc_) << "\n"
  //      << "max(dd) = "<< (dd.weightedAverage(mu.mesh().Vsc())).value() << "\n"
  //      << "max(rhoc/rho()) = "<< max(rhoc_/rho()).value() << "\n"
  //      << "max(rho) = "<< max(rho()).value() << " min(rho) = " << min(rho()).value() << "\n"
  //      << "min(mu) = "<< min(mu).value() << " max(mu) = "<< max(mu).value() << endl;


}


// ************************************************************************* //
