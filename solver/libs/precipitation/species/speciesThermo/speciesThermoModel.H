/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::speciesThermoModel

Description
    Base class for speciesThermo classes that encompass speciation of a set of
    reacting species.

SourceFiles
    speciesThermoModel.C

\*---------------------------------------------------------------------------*/

#ifndef speciesThermoModel_H
#define speciesThermoModel_H

#include "primitiveFieldsFwd.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "multivariateScheme.H"

#include "speciesMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class speciesThermoModel Declaration
\*---------------------------------------------------------------------------*/

class speciesThermoModel
{

protected:

    // Protected data
  
        const fvMesh& mesh_;  
   
        volScalarField& rho_;

        const dictionary& speciesDict_;

        speciesMixture& speciesMixture_;

        Switch printParams_;

    // Private Member Functions

        //- Disallow copy construct
        speciesThermoModel(const speciesThermoModel&);

        //- Disallow default bitwise assignment
        void operator=(const speciesThermoModel&);


public:

    //- Runtime type information
    TypeName("speciesThermoModel");


    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            speciesThermoModel,
            speciesThermoModel,
            (
                volVectorField& U,
                volScalarField& rho,
                dictionary& speciesDict
            ),
            (U, rho, speciesDict)
        );


    // Constructors

        //- Construct from components
        speciesThermoModel
        ( 
            volVectorField& U,
            volScalarField& rho,
            dictionary& speciesDict
        );

    // Selectors
   
        static autoPtr<speciesThermoModel> New
        (
            volVectorField& U,
            volScalarField& rho,
            dictionary& speciesDict
        );  

    //- Destructor
    virtual ~speciesThermoModel();


    // Member Functions

        // Evolution
            
            //- Source term in species equation
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;        

            //- Return ionic strength
            virtual volScalarField I() const = 0;

            //- Return ionic activity
            virtual volScalarField gamma(const word& specie) const = 0;

            //- Return Debye-Huckel constant
            //virtual scalar DebyeHuckel() const = 0;
            
            //- Ionic activity product
            virtual volScalarField IAP() const = 0;

            //- Supersaturation ratio
            virtual volScalarField SR() const = 0;

            //- Saturation index
            virtual volScalarField SI() const = 0;
            
            //- Absolute saturation
            virtual volScalarField S() const = 0;

            //- Relative saturation 
            virtual volScalarField relS() const = 0;

            //- Activity based supersaturation
            virtual volScalarField Sa() const = 0;



        // Lookup functions
            
            //- Return mesh
            inline const fvMesh& mesh() const;
            
            //- Return mixture density
            inline const volScalarField& rho() const;
            inline volScalarField& rho();

            //- speciesDict getters
            inline const dictionary& speciesDict() const;
          //  inline dictionary& speciesDict();

            //- speciesMixture getters
            inline const speciesMixture& speciesComposition() const;
            inline speciesMixture& speciesComposition();


        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read() = 0;

        //- Print concentrations of reactants
        virtual void printThermoParams() = 0;
        
        virtual void updateConcentrations() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
#include "speciesThermoModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
