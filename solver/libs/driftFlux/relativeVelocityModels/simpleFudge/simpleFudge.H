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
    Foam::simpleFudge

Description
    Simple relative velocity model

SourceFiles
    simpleFudge.C

\*---------------------------------------------------------------------------*/

#ifndef simpleFudge_H
#define simpleFudge_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class simpleFudge Declaration
\*---------------------------------------------------------------------------*/

class simpleFudge
:
    public relativeVelocityModel
{
    // Private data

        //- A coefficient
        dimensionedScalar a_;

        //- Drift velocity
        dimensionedVector V0_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;

        //- Terminal velocity field
        volVectorField Vt_;  //TODO FUDGE

        //- Gravity const vector
        dimensionedVector g_; //TODO FUDGE 


public:

    //- Runtime type information
    TypeName("simpleFudge");


    // Constructors

        //- Construct from components
        simpleFudge
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~simpleFudge();


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
