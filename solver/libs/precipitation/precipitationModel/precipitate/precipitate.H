/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::incompressible::precipitate

Description
    Turbulence model for precipitate incompressible flow.

SourceFiles
    precipitate.C

\*---------------------------------------------------------------------------*/

#ifndef precipitate_H
#define precipitate_H

#include "precipitationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class precipitate Declaration
\*---------------------------------------------------------------------------*/

class precipitate
:
    public precipitationModel
{

public:

    //- Runtime type information
    TypeName("precipitate");

    // Constructors

        //- Construct from components
        precipitate
        (
            volVectorField& U,
            volScalarField& rho,
            surfaceScalarField& phi
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<precipitate> New
        (
            volVectorField& U,
            volScalarField& rho,
            surfaceScalarField& phi
        );


    //- Destructor
    virtual ~precipitate()
    {}


    // Member Functions
        
        // Evolution
        
            //- Species source term 
            virtual tmp<volScalarField> precipitationSource(volScalarField& Y) const;
          //  virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

          //  virtual volScalarField& growthSource();

         //   virtual void solvePBE();

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

