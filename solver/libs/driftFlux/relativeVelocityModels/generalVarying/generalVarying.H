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

Class
    Foam::generalVarying

Description
    General relative velocity model

SourceFiles
    generalVarying.C

\*---------------------------------------------------------------------------*/

#ifndef generalVarying_H
#define generalVarying_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class generalVarying Declaration
\*---------------------------------------------------------------------------*/

class generalVarying
:
    public relativeVelocityModel
{
    // Private data

        //- A coefficient
        dimensionedScalar a_;

        //- a1 coefficient
        dimensionedScalar a1_;

        //- Drift velocity
        dimensionedVector V0_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;

        //- Drift velocity field
        volVectorField Vt_;

        dimensionedVector g_; //TODO FUDGE


public:

    //- Runtime type information
    TypeName("generalVarying");


    // Constructors

        //- Construct from components
        generalVarying
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~generalVarying();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
