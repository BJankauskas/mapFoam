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
    Foam::precipitationModel

Description
    Base-class for all transport models used by the incompressible turbulence
    models.

SourceFiles
    precipitationModel.C

\*---------------------------------------------------------------------------*/

#ifndef precipitationModel_H
#define precipitationModel_H

#include "primitiveFieldsFwd.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "multivariateScheme.H"

#include "speciesThermoModel.H"

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class precipitationModel Declaration
\*---------------------------------------------------------------------------*/

class precipitationModel
:
    public IOdictionary 
{

protected:

    // Protected data
  
        const fvMesh& mesh_;  

      //  const word precipitationModelName_;
      //  const word alphaPrecipitationModelName_;

        dictionary& speciesDict_;
   
        autoPtr<speciesThermoModel> speciesThermoModel_;

        IOdictionary populationBalanceProperties_;

        autoPtr<populationBalanceModel> populationBalanceModel_;

        Switch addAlphaSource_;

        const word precipitant_;

    // Private Member Functions

        //- Disallow copy construct
        precipitationModel(const precipitationModel&);

        //- Disallow default bitwise assignment
        void operator=(const precipitationModel&);


public:

    //- Runtime type information
    TypeName("precipitationModel");


    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            precipitationModel,
            precipitationModel,
            (
                volVectorField& U,
                volScalarField& rho,
                surfaceScalarField& phi
            ),
            (U, rho, phi)
        );


    // Constructors

        //- Construct from components
        precipitationModel
        ( 
            volVectorField& U,
            volScalarField& rho,
            surfaceScalarField& phi
        );

    // Selectors
   
        static autoPtr<precipitationModel> New
        (
            volVectorField& U,
            volScalarField& rho,
            surfaceScalarField& phi
        );  

    //- Destructor
    virtual ~precipitationModel();


    // Member Functions

        // Evolution
            
            //- Source term in species equation
            virtual tmp<volScalarField> precipitationSource(volScalarField& Y) const = 0;
            virtual tmp<volScalarField> alphaPrecipitationSource() const = 0;
            //virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;        

        // Lookup functions
            
            //- Return mesh
            inline const fvMesh& mesh() const;

            inline const word& precipitationModelName() const;
            inline const word& alphaPrecipitationModelName() const;

            //- speciesDict getters
            inline const dictionary& speciesDict() const;
            //inline dictionary& speciesDict();

            //- speciesMixture getters
            inline const speciesThermoModel& speciesThermo() const;
            inline speciesThermoModel& speciesThermo();

            //- populationBalanceModel getters
            inline const IOdictionary& populationBalanceProperties() const;
            inline const populationBalanceModel& populationBalance() const;
      //      virtual populationBalanceModel& populationBalance():
      
            inline const word& precipitant() const; 

            //virtual volScalarField& growthSource() = 0;

            void solvePBE();


        //- Read transportProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "precipitationModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

