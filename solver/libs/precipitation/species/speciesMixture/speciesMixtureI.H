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
    Foam::speciesMixture

SourceFiles
    speciesMixture.C

\*---------------------------------------------------------------------------*/

inline const Foam::speciesTable&
Foam::speciesMixture::species() const
{
    return species_;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::speciesMixture::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::speciesMixture::Y() const
{
    return Y_;
}


inline Foam::volScalarField& Foam::speciesMixture::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::speciesMixture::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::speciesMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::speciesMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}


inline bool Foam::speciesMixture::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}

/******* Molecular weight getters ********/

inline Foam::PtrList<Foam::dimensionedScalar>&
Foam::speciesMixture::W()
{
    return W_;
}

inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::speciesMixture::W() const 
{
    return W_;
}


inline Foam::dimensionedScalar&
Foam::speciesMixture::W
(
    const label i
)
{
    return W_[i];
}


inline const Foam::dimensionedScalar&
Foam::speciesMixture::W
(
    const label i
) const
{
    return W_[i];
}

inline Foam::dimensionedScalar&
Foam::speciesMixture::W
(
    const word& specieName
)
{
    return W_[species_[specieName]];
}

inline const Foam::dimensionedScalar&
Foam::speciesMixture::W
(
    const word& specieName
) const
{
    return W_[species_[specieName]];
}

/****** IONIC CHARGES getters  ****/

inline Foam::List<Foam::scalar>&
Foam::speciesMixture::Z()
{
    return Z_;
}

inline const Foam::List<Foam::scalar>&
Foam::speciesMixture::Z() const 
{
    return Z_;
}


inline Foam::scalar&
Foam::speciesMixture::Z
(
    const label i
)
{
    return Z_[i];
}


inline const Foam::scalar&
Foam::speciesMixture::Z
(
    const label i
) const
{
    return Z_[i];
}

inline Foam::scalar&
Foam::speciesMixture::Z
(
    const word& specieName
)
{
    return Z_[species_[specieName]];
}

inline const Foam::scalar&
Foam::speciesMixture::Z
(
    const word& specieName
) const
{
    return Z_[species_[specieName]];
}

// ************************************************************************* //
