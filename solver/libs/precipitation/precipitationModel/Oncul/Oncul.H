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
    Foam::incompressible::Oncul

Description
    Turbulence model for Oncul incompressible flow.

SourceFiles
    Oncul.C

\*---------------------------------------------------------------------------*/

#ifndef Oncul_H
#define Oncul_H

#include "precipitationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Oncul Declaration
\*---------------------------------------------------------------------------*/

class Oncul
:
    public precipitationModel
{
private:

    // Growth rate and exponent
    dimensionedScalar Cg_;
    scalar Ng_;

    // Shape factor
    scalar kv_;
    
    // Precipitant name
    word precipitant_;

    // Product species molar weight
    dimensionedScalar Ws_;

    // Crystal density
    const volScalarField& rhod_;

public:

    //- Runtime type information
    TypeName("Oncul");

    // Constructors

        //- Construct from components
        Oncul
        (
            volVectorField& U,
            volScalarField& rho,
            surfaceScalarField& phi
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<Oncul> New
        (
            volVectorField& U,
            volScalarField& rho,
            surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Oncul()
    {}


    // Member Functions
        
        // Evolution
        
            //- Species source term 
            virtual tmp<volScalarField> precipitationSource(volScalarField& Y) const;
            virtual tmp<volScalarField> alphaPrecipitationSource() const;
          //  virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

          //  virtual volScalarField& growthSource();

           // virtual void solvePBE();

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

