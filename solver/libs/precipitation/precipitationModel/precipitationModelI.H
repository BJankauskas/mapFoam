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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::precipitationModel::mesh() const
{
    return mesh_;
}

//inline const Foam::word& Foam::precipitationModel::precipitationModelName() const
//{
//    return precipitationModelName_;
//}

//inline const Foam::word& Foam::precipitationModel::alphaPrecipitationModelName() const
//{
//    return alphaPrecipitationModelName_;
//}

inline const Foam::dictionary& Foam::precipitationModel::speciesDict() const
{
    return speciesDict_;
}

//inline Foam::dictionary& Foam::precipitationModel::speciesDict() 
//{
//    return speciesDict_;
//}

inline const Foam::speciesThermoModel& 
Foam::precipitationModel::speciesThermo() const
{
    return speciesThermoModel_();
}

inline Foam::speciesThermoModel& Foam::precipitationModel::speciesThermo() 
{
    return speciesThermoModel_();
}

inline const Foam::IOdictionary&
Foam::precipitationModel::populationBalanceProperties() const
{
    return populationBalanceProperties_;
};

inline const Foam::populationBalanceModel& 
Foam::precipitationModel::populationBalance() const
{
    return populationBalanceModel_();
}

inline const Foam::word&
Foam::precipitationModel::precipitant() const
{
    return precipitant_;
}

/*
inline const Foam::volScalarField& speciesSource() const
{
    return speciesSource_;
}
*/
// ************************************************************************* //
