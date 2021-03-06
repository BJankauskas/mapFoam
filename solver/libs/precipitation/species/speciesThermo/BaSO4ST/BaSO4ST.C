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

\*---------------------------------------------------------------------------*/

#include "BaSO4ST.H"
#include "Time.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BaSO4ST, 0);
addToRunTimeSelectionTable(speciesThermoModel, BaSO4ST, speciesThermoModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BaSO4ST::BaSO4ST
(
    volVectorField& U,
    volScalarField& rho,
    dictionary& speciesDict
)
:
    speciesThermoModel(U,rho, speciesDict),
    coeffsDict_(this->speciesDict().subDict(this->type() + "Coeffs")),
    pH_(readScalar(coeffsDict_.lookup("pH"))),
    gamma_(readScalar(coeffsDict_.lookup("gamma"))),
   // kv_(readScalar(coeffsDict_.lookup("kv"))),
   // ka_(readScalar(coeffsDict_.lookup("ka"))),

    // Equilibrium constants
    equilibriumCoeffs_(coeffsDict_.subDict("equilibriumConstants")),
    Ksp_
    (
        "Ksp",
        dimensionSet(0,-6,0,0,2,0,0), 
        readScalar(equilibriumCoeffs_.lookup("Ksp"))
    ),

    // Concentration fields
    
    C_(this->speciesComposition().species().size()),
    Ct_(this->speciesComposition().species().size())
    
{
    // Populate C_ and Ct_
    
    const speciesTable& species_ = this->speciesComposition().species();

    forAll(species_, i)
    {
        // NOTE Ommit header check as used in species mixture class
        C_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "C_" + species_[i],
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
            )
        );   

        Ct_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Ct_" + species_[i],
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
            )
        );       
    }

    // Initialise concentrations fields
    this->updateConcentrations();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<BaSO4ST> BaSO4ST::New
(
    volVectorField& U,
    volScalarField& rho,
    dictionary& speciesDict
)
{
    return autoPtr<BaSO4ST>
    (
        new BaSO4ST(U, rho, speciesDict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::BaSO4ST::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimMass/dimTime)
    );

    return tSu;
} 

Foam::dimensionedScalar Foam::BaSO4ST::Hconc() const
{
    dimensionedScalar Hc_
    (
        "Hc_",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,-pH_)   
    );
        
    return Hc_;
}

Foam::dimensionedScalar Foam::BaSO4ST::OHconc() const
{
    dimensionedScalar OHc_
    (   
        "OHc_",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10, -(14-pH_))
    );

    return OHc_;
}       

/****************************************************/
/*               Ion concentrations                 */
/****************************************************/

Foam::volScalarField Foam::BaSO4ST::Ct(const word& specie) const
{
    // Mass fractions
    const volScalarField& Yi = this->speciesComposition().Y(specie);

    // Mol weights
    const dimensionedScalar& Wi = this->speciesComposition().W(specie);
    
    // Mixture density
    const volScalarField& rho = this->rho();
     
    volScalarField Ct
    {
        rho*Yi/Wi 
    };

    return Ct;
}

//Foam::volScalarField Foam::BaSO4ST::C_Ba() const
//{
//    // Mass fractions (C_Ts in Ye et al 2018)
//    const volScalarField& Y_Ba = this->speciesComposition().Y("Ba");
//    
//    // Mol Weights
//    const dimensionedScalar& W_Ba = this->speciesComposition().W("Ba");
//
//    // Mixture density
//    const volScalarField& rho = this->rho();
//   
//    volScalarField tC_Ba
//    (
//        rho*Y_Ba/W_Ba
//    );
//
//    return tC_Ba;
//}
//
//Foam::volScalarField Foam::BaSO4ST::C_SO4() const
//{
//    // Mass fraction
//    const volScalarField& Y_SO4= this->speciesComposition().Y("SO4");
//
//    // Mol weight
//    const dimensionedScalar& W_SO4= this->speciesComposition().W("SO4");
//
//    // Mixture density
//    const volScalarField& rho = this->rho();
//
//    volScalarField tC_SO4
//    (
//        (rho*Y_SO4/W_SO4)
//    );
//
//    return tC_SO4;
//}
//
//Foam::volScalarField Foam::BaSO4ST::C_BaSO4() const
//{
//    // Mass fraction
//    const volScalarField& Y_BaSO4= this->speciesComposition().Y("BaSO4");
//
//    // Mol weight
//    const dimensionedScalar& W_BaSO4= this->speciesComposition().W("BaSO4");
//
//    // Mixture density
//    const volScalarField& rho = this->rho();
//
//    volScalarField tC_BaSO4
//    (
//        (rho*Y_BaSO4/W_BaSO4)
//    );
//
//    return tC_BaSO4;
//}

/****************************************************/
/*             Ionic activity functions             */
/****************************************************/

Foam::volScalarField
Foam::BaSO4ST::I() const
{

   volScalarField tI
   (
       IOobject
       (
           "I", 
           this->mesh().time().timeName(),
           this->mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE,
           false
       ),
       this->mesh(),
       dimensionedScalar("zero",dimMoles/dimVolume,0)
   );


   // TODO In the future use PtrList for Ct fields

//   forAll(C_, i)
//   {
//        const scalar& Zi = this->speciesComposition().Z(i);
//        volScalarField& Ci = C_[i];
//
//        I += pow(Zi,2) * Ci;
//   }
//
//   I.dimensions().reset(dimless);
    
   return tI;

}

Foam::volScalarField
Foam::BaSO4ST::gamma(const word& specie) const
{

    volScalarField tGamma
    (
         IOobject
         (
             "gamma", 
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false
         ),
         this->mesh(),
         dimensionedScalar("zero", dimless, gamma_)
    );
 
//    volScalarField& Gamma = tGamma.ref();
//
//    scalar ADH = readScalar(this->speciesDict().lookup("ADH"));
//    scalar Zi = speciesThermo().Z(Y.name());
//
//    volScalarField& I = this->I();
//
//    Gamma = ADH*pow(Zi,2)*(pow(I, 0.5)/(1 + pow(I, 0.5))) - 0.3 * I; 
//    
//    return exp(-Gamma);

      return tGamma;
}




/****************************************************/
/*               Saturation Indexes                 */
/****************************************************/

Foam::volScalarField Foam::BaSO4ST::IAP() const
{
    return volScalarField::null();
}

// Supersaturation ratio
Foam::volScalarField Foam::BaSO4ST::SR() const
{
    volScalarField tSR
    (
        IOobject
        (
            "tSR",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar("unit", dimless, 1.0)
    );

    const speciesTable& species = this->speciesComposition().species();

    forAll(species, i)
    {
        const scalar rExp = 
                readScalar(this->speciesDict().subDict(species[i]).lookup("rExponent"));

        volScalarField Ci(this->Ct(species[i]));

        tSR *=  pow(Ci, rExp);

    }

    // Info<<"min(IAP) = " << min(tSR).value() << endl;
    // Info<<"1/Ksp = " << (1/Ksp_).value() << endl;

    tSR /= Ksp_;

    tSR.max(0);

    return tSR;
}

// Supersaturation index
Foam::volScalarField Foam::BaSO4ST::SI() const
{
    volScalarField tSI(this->SR());
    tSI.max(1.0);
    tSI = log10(tSI);

    return tSI;
}

// Absolute supersaturation 
Foam::volScalarField Foam::BaSO4ST::S() const
{
    volScalarField tS(this->SR());

    tS = pow(Ksp_.value(), 0.5) * (pow(tS, 0.5) - 1);
    tS.max(0);
   
    // tS.dimensions().reset(dimensionSet(0,-3,0,0,1,0,0));

    return tS;
}

// Relative supersaturation
Foam::volScalarField Foam::BaSO4ST::relS() const
{
    volScalarField tRelS(this->SR());

    tRelS = pow(tRelS, 0.5) - 1;
    tRelS.max(0);

    return tRelS;
}

// Activity supersaturation or sqrt(Supersaturation ratio)
Foam::volScalarField Foam::BaSO4ST::Sa() const
{
    volScalarField tSa(this->SR());

    tSa = gamma_ * pow(tSa, 0.5);

    return tSa;
}

void BaSO4ST::correct()
{}


bool BaSO4ST::read()
{
    return true;
}
void BaSO4ST::printThermoParams()
{
    if(printParams_)
    {
        Info<<"/***************************************************************/ \n" 
        <<"/**********   Supersaturation parameters output: ***************/ \n"
        <<"/***************************************************************/ \n"

        <<"Saturation ratio min(SR) = "<< min(this->SR()).value()
        <<" max(SR) = "<< max(this->SR()).value() << "\n"

        <<"Saturation Index min(SI) = "<< min(this->SI()).value()
        <<" max(SI) = "<< max(this->SI()).value() << "\n"

        <<"Absolute saturation min(S) = "<< min(this->S()).value()
        <<" max(S) = "<< max(this->S()).value() << "\n"

        <<"cbrt(SR) min(Sa) = "<< min(this->Sa()).value()
        <<" max(Sa) = "<< max(this->Sa()).value() <<  "\n"

        <<"Relative supersaturation min(relS) = "<< min(this->relS()).value()
        <<" max(relS) = "<< max(this->relS()).value() << endl;

        scalar reqH2O = 0;

        forAll(this->C(), i)
        {   
            volScalarField& Ci = this->C(i);
            dimensionedScalar& Wi = this->speciesComposition().W(i);

            const word& specie = this->speciesComposition().Y(i).name();
            const scalar reqC = this->speciesDict().subDict(specie).lookupOrDefault<scalar>("requiredC", 0.0);

            if(specie != "H2O")
            {
                Info<<"min("<< specie << ") = "<< min(Ci).value() 
                    <<" max("<< specie << ") = "<< max(Ci).value() << "\n"
                    <<"Required(Y_"<< specie << ") = " << reqC*Wi.value()/max(this->rho()).value() << endl; 
                    reqH2O += reqC*Wi.value()/max(this->rho()).value(); 
            } else {
            
                Info<<"min("<< specie << ") = "<< min(Ci).value() 
                    <<" max("<< specie << ") = "<< max(Ci).value() << "\n"
                    <<"Required(Y_"<< specie << ") = " << 1 - reqH2O << endl;        
            }
        }

      //  Info<<"Chemical Conversion = "<< max(C("MAP")).value()/(0.04 - pow(Ksp_.value(), 0.33)) << endl;

        Info<<"/*************************************************************/" << endl;
    }
}

void BaSO4ST::updateConcentrations()
{
    forAll(this->C(), i)
    {  
        // Update total concentration fields from mass fractions
        volScalarField& Ci = this->C(i);
        volScalarField& Yi = this->speciesComposition().Y(i);
        dimensionedScalar& Wi = this->speciesComposition().W(i);
        volScalarField& rho = this->rho();

        Ci = rho * Yi / Wi;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
