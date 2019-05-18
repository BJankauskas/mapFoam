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

Class
    Foam::populationBalanceSubModels::nucleationModels::powerLawNucleation

Description
    powerLawNucleation nucleation model for soot formation.

SourceFiles
    powerLawNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawNucleation_H
#define powerLawNucleation_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                        Class powerLawNucleation Declaration
\*---------------------------------------------------------------------------*/

class powerLawNucleation
:
    public nucleationModel
{
private:

    // Private data
    
        dictionary coeffsDict_;

        //- Nucleation rate
        dimensionedScalar Kg_;

        //- Nucleation activation parameter
        scalar exponent_;

protected:

        volScalarField J_;

public:

        //- Runtime type information
        TypeName("powerLawNucleation");


    // Constructors

        //- Construct from components
        powerLawNucleation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~powerLawNucleation();


    // Member Functions
    
        //- Induction time
        void tInduction
        (
            const volUnivariateMoment& moment
        ) const;

        //- Nucleation source
        virtual tmp<volScalarField> nucleationSource
        (
            const volUnivariateMoment& moment
        ); //const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
