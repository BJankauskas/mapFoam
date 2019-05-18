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
    Foam::MAPST

Description

SourceFiles
    MAPST.C

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

/*          Access to total concentrations      */

inline Foam::PtrList<Foam::volScalarField>&
Foam::MAPST::Ct()
{
    return Ct_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::MAPST::Ct() const
{
    return Ct_;
}


inline Foam::volScalarField& Foam::MAPST::Ct(const label i)
{
    return Ct_[i];
}


inline const Foam::volScalarField& Foam::MAPST::Ct
(
    const label i
) const
{
    return Ct_[i];
}


inline Foam::volScalarField& Foam::MAPST::Ct
(
    const word& specieName
)
{
    const speciesTable& species = this->speciesComposition().species();
    return Ct_[species[specieName]];
}


inline const Foam::volScalarField& Foam::MAPST::Ct
(
    const word& specieName
) const
{
    const speciesTable& species = this->speciesComposition().species();
    return Ct_[species[specieName]];
}

/*       Access to ion concentrations           */

inline Foam::PtrList<Foam::volScalarField>&
Foam::MAPST::C()
{
    return C_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::MAPST::C() const
{
    return C_;
}


inline Foam::volScalarField& Foam::MAPST::C(const label i)
{
    return C_[i];
}


inline const Foam::volScalarField& Foam::MAPST::C
(
    const label i
) const
{
    return C_[i];
}


inline Foam::volScalarField& Foam::MAPST::C
(
    const word& specieName
)
{
    const speciesTable& species = this->speciesComposition().species();
    return C_[species[specieName]];
}


inline const Foam::volScalarField& Foam::MAPST::C
(
    const word& specieName
) const
{
    const speciesTable& species = this->speciesComposition().species();
    return C_[species[specieName]];
}

// ************************************************************************* //
