/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::speciesMixture

Description
    Class containing information on all reacting species. 
    Based on Foam::basicMultiComponentMixture.

SourceFiles
    speciesMixture.C

\*---------------------------------------------------------------------------*/

#ifndef speciesMixture_H
#define speciesMixture_H

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"
#include "typeInfo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class speciesMixture Declaration
\*---------------------------------------------------------------------------*/

class speciesMixture
{

    // Private Member Functions

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Construct as copy (not implemented)
        speciesMixture(const speciesMixture&);

protected:

    // Protected data
   
        //- Table of specie names
        speciesTable species_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;

        //- Species molecular weights
        PtrList<dimensionedScalar> W_;

        //- Species ionic charge
        List<scalar> Z_;


public:

    // Constructors

        //- Construct from dictionary, specie names and  mesh //TODO might need to include phaseName for later on
        speciesMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&
        );

        //- Construct from dictionary and  mesh //TODO might need to include phaseName for later on
        speciesMixture(const dictionary&, const fvMesh&);

    //- Destructor
    virtual ~speciesMixture()
    {}


    // Member functions

        // Access
       
            //- Return the table of species
            inline const speciesTable& species() const;

            //- Return the mass-fraction fields
            inline PtrList<volScalarField>& Y();

            //- Return the const mass-fraction fields
            inline const PtrList<volScalarField>& Y() const;

            //- Return the mass-fraction field for a specie given by index
            inline volScalarField& Y(const label i);

            //- Return the const mass-fraction field for a specie given by index
            inline const volScalarField& Y(const label i) const;

            //- Return the mass-fraction field for a specie given by name
            inline volScalarField& Y(const word& specieName);

            //- Return the const mass-fraction field for a specie given by name
            inline const volScalarField& Y(const word& specieName) const;

            //- Does the mixture include this specie?
            inline bool contains(const word& specieName) const;

            //- Return the molecular weights
            inline PtrList<dimensionedScalar>& W();

            //- Return the const molecular weights
            inline const PtrList<dimensionedScalar>& W() const;

            //- Return the molecular weights for a specie given by index
            inline dimensionedScalar& W(const label i);

            //- Return the const molecular weights for a specie given by index
            inline const dimensionedScalar& W(const label i) const;

            //- Return the molecular weights for a specie given by name
            inline dimensionedScalar& W(const word& specieName);

            //- Return the const molecular weights for a specie given by name
            inline const dimensionedScalar& W(const word& specieName) const;

            //- Return the ionic charges
            inline List<scalar>& Z();

            //- Return the const ionic charges
            inline const List<scalar>& Z() const;

            //- Return the ionic charge for a specie given by index
            inline scalar& Z(const label i);

            //- Return the const ionic charge for a specie given by index
            inline const scalar& Z(const label i) const;

            //- Return the ionic charge for a specie given by name
            inline scalar& Z(const word& specieName);

            //- Return the const ionic charge for a specie given by name
            inline const scalar& Z(const word& specieName) const;


        //- Read dictionary
        void read(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "speciesMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
