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
    Foam::incompressible::MAPST

Description
    Turbulence model for MAPST incompressible flow.

SourceFiles
    MAPST.C

\*---------------------------------------------------------------------------*/

#ifndef MAPST_H
#define MAPST_H

#include "speciesThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MAPST Declaration
\*---------------------------------------------------------------------------*/

class MAPST
:
    public speciesThermoModel
{
protected:

    // Protected data
        
        const dictionary& coeffsDict_;

        // pH level
        scalar pH_;

        bool gammaVariation_;
        scalar gamma_;
        scalar ADH_;

        // Equilibrium constants
        const dictionary& equilibriumCoeffs_;

        dimensionedScalar kMgOH_;
        dimensionedScalar kMgHPO4_;
        dimensionedScalar kMgH2PO4_;
        dimensionedScalar kMgPO4_;
        dimensionedScalar kH3PO4_;
        dimensionedScalar kH2PO4_;
        dimensionedScalar kHPO4_;
        dimensionedScalar kNH4_;
        dimensionedScalar kH2O_;

        dimensionedScalar Ksp_;

        // Concentration fields
        PtrList<volScalarField> C_;
        PtrList<volScalarField> Ct_;
        
        /*
        volScalarField C_PO4_;
        volScalarField C_NH4_;
        volScalarField C_MG_;
        volScalarField C_MAP_;
        */

/*
        volScalarField Ct_PO4_;
        volScalarField Ct_NH4_;
        volScalarField Ct_MG_;
        volScalarField Ct_MAP_;
        */
    

public:

    //- Runtime type information
    TypeName("MAPST");

    // Constructors

        //- Construct from components
        MAPST
        (
            volVectorField& U,
            volScalarField& rho,
            dictionary& speciesDict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<MAPST> New
        (
            volVectorField& U,
            volScalarField& rho,
            dictionary& speciesDict
        );


    //- Destructor
    virtual ~MAPST()
    {}


    // Member Functions
        
        // Evolution
        
            //- Species source term 
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            dimensionedScalar Hconc() const;
            dimensionedScalar OHconc() const;

            //- Speciation coefficients
            volScalarField K0() const;
            volScalarField K1() const;
            volScalarField K2() const;
            volScalarField K3() const;
    
            //- Calculate ion concentrations
            volScalarField C_PO4() const;
            volScalarField C_NH4() const;
            volScalarField C_MG() const;
            volScalarField C_MAP() const;

            //- Calculate total concentrations from mass-fractions
       //     volScalarField Ct_PO4() const;
       //     volScalarField Ct_NH4() const;
       //     volScalarField Ct_MG() const;
       //     volScalarField Ct_MAP() const;

            //- Return ionic strength
            virtual volScalarField I() const;

            //- Return ionic activity
            virtual volScalarField gamma(const word& specie) const;

            //- Ionic activity product
            virtual volScalarField IAP() const;    

            //- Supersaturation ratio
            virtual volScalarField SR() const;

            //- Saturation index
            virtual volScalarField SI() const;
            
            //- Absolute saturation
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
            
        //- Correct the MAPST viscosity
        virtual void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read();

        virtual void printThermoParams();
        virtual void updateConcentrations();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MAPSTI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
