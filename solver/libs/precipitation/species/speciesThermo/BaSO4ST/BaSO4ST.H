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
    Foam::incompressible::BaSO4ST

Description
    Turbulence model for BaSO4ST incompressible flow.

SourceFiles
    BaSO4ST.C

\*---------------------------------------------------------------------------*/

#ifndef BaSO4ST_H
#define BaSO4ST_H

#include "speciesThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BaSO4ST Declaration
\*---------------------------------------------------------------------------*/

class BaSO4ST
:
    public speciesThermoModel
{
protected:

    // Protected data
        
        const dictionary& coeffsDict_;

        // pH level
        scalar pH_;
        scalar gamma_;

        // Shape factors
        //scalar kv_;
        //scalar ka_;

        // Equilibrium constants
        const dictionary& equilibriumCoeffs_;

        dimensionedScalar Ksp_;

        volScalarField C_SO4_;
        volScalarField C_Ba_;
        volScalarField C_BaSO4_;
    

public:

    //- Runtime type information
    TypeName("BaSO4ST");

    // Constructors

        //- Construct from components
        BaSO4ST
        (
            volVectorField& U,
            volScalarField& rho,
            dictionary& speciesDict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<BaSO4ST> New
        (
            volVectorField& U,
            volScalarField& rho,
            dictionary& speciesDict
        );


    //- Destructor
    virtual ~BaSO4ST()
    {}


    // Member Functions
        
        // Evolution
        
            //- Species source term 
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            dimensionedScalar Hconc() const;
            dimensionedScalar OHconc() const;


            //- Calculate ion concentrations
            //volScalarField C_Ba() const;
            //volScalarField C_SO4() const;
            //volScalarField C_BaSO4() const;
            
            volScalarField Ct(const word& specie) const;

            //- Return ionic strength
            virtual volScalarField I() const;
           
            //- Return ionic activity coefficient
            virtual volScalarField gamma(const word& specie) const;            

            //- Ionic activity product
            virtual volScalarField IAP() const;    

            //- Supersaturation ratio
            virtual volScalarField SR() const;

            //- Saturation index
            virtual volScalarField SI() const;

            //-Absolute saturation
            virtual volScalarField S() const;
            
            //- Relative saturation
            virtual volScalarField relS() const;

            //- Activity based supersaturation
            virtual volScalarField Sa() const;
            
        //- Access functions
            
            //- Return the total concentration fields
            inline PtrList<volScalarField>& Ct();

            //- Return the const total concentration fields
            inline const PtrList<volScalarField>& Ct() const;

            //- Return the total concentration field for a specie given by index
            inline volScalarField& Ct(const label i);

            //- Return the const total concentration field for a specie given by index
            inline const volScalarField& Ct(const label i) const;

            //- Return the total concentration field for a specie given by name
            inline volScalarField& Ct(const word& specieName);

            //- Return the const total concentration field for a specie given by name
            inline const volScalarField& Ct(const word& specieName) const;

            //- Return the ion concentration fields
            inline PtrList<volScalarField>& C();

            //- Return the const ion concentration fields
            inline const PtrList<volScalarField>& C() const;

            //- Return the ion concentration field for a specie given by index
            inline volScalarField& C(const label i);

            //- Return the const ion concentration field for a specie given by index
            inline const volScalarField& C(const label i) const;

            //- Return the ion concentration field for a specie given by name
            inline volScalarField& C(const word& specieName);

            //- Return the const ion concentration field for a specie given by name
            inline const volScalarField& C(const word& specieName) const;        


        //- Correct the BaSO4ST viscosity
        virtual void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read();

        virtual void printThermoParams();
        virtual void updateConcentrations();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
