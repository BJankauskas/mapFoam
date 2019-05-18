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
    Base-class for all transport models used by the incompressible turbulence
    models.

SourceFiles
    speciesThermoModel.C

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::speciesThermoModel::mesh() const
{
    return mesh_;
}

inline const Foam::volScalarField& Foam::speciesThermoModel::rho() const
{
    return rho_;
}

inline Foam::volScalarField& Foam::speciesThermoModel::rho()
{
    return rho_;
}

inline Foam::speciesMixture& Foam::speciesThermoModel::speciesComposition()
{
    return speciesMixture_;
}

inline const Foam::speciesMixture& Foam::speciesThermoModel::speciesComposition() const
{
    return speciesMixture_;
}

//inline Foam::dictionary& Foam::speciesThermoModel::speciesDict()
//{
//    return speciesDict_;
//}

inline const Foam::dictionary& Foam::speciesThermoModel::speciesDict() const
{
    return speciesDict_;
}

// ************************************************************************* //
