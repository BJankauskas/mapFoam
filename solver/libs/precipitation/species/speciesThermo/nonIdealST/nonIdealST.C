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

#include "nonIdealST.H"
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

defineTypeNameAndDebug(nonIdealST, 0);
addToRunTimeSelectionTable(speciesThermoModel, nonIdealST, speciesThermoModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonIdealST::nonIdealST
(
    volVectorField& U,
    volScalarField& rho,
    dictionary& speciesDict
)
:
    speciesThermoModel(U,rho, speciesDict),
    coeffsDict_(this->speciesDict().subDict(this->type() + "Coeffs")),
    pH_(readScalar(coeffsDict_.lookup("pH"))),
    gammaVariation_(coeffsDict_.lookupOrDefault<Switch>("gammaVariation", false)),
    gamma_(readScalar(coeffsDict_.lookup("gamma"))),
    ADH_(coeffsDict_.lookupOrDefault<scalar>("ADH", 1)),
   // kv_(readScalar(coeffsDict_.lookup("kv"))),
   // ka_(readScalar(coeffsDict_.lookup("ka"))),

    // Equilibrium constants
    equilibriumCoeffs_(coeffsDict_.subDict("equilibriumConstants")),
    nReactants_(readScalar(coeffsDict_.lookup("nReactants"))),
    Ksp_
    (
        "Ksp",
        pow((dimMoles/dimVolume), nReactants_),
        readScalar(equilibriumCoeffs_.lookup("Ksp"))
        //pow(10, readScalar(equilibriumCoeffs_.lookup("pKsp")))
    ),
    
    // Ion concentration fields
    
    C_(this->speciesComposition().species().size())
{

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
    }

    // Initialise concentration fields
    this->updateConcentrations();

}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<nonIdealST> nonIdealST::New
(
    volVectorField& U,
    volScalarField& rho,
    dictionary& speciesDict
)
{
    return autoPtr<nonIdealST>
    (
        new nonIdealST(U, rho, speciesDict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::nonIdealST::R
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

Foam::dimensionedScalar Foam::nonIdealST::Hconc() const
{
    dimensionedScalar Hc_
    (
        "Hc_",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,-pH_)   
    );
        
    return Hc_;
}

Foam::dimensionedScalar Foam::nonIdealST::OHconc() const
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

Foam::volScalarField Foam::nonIdealST::Ct(const word& specie) const
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

/****************************************************/
/*             Ionic activity functions             */
/****************************************************/

Foam::volScalarField
Foam::nonIdealST::I() const
{

   volScalarField I
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

   const speciesTable& species = speciesComposition().species(); 

   forAll(species, i)
   {
        const scalar& Zi = this->speciesComposition().Z(i);
        const volScalarField& Ci = this->C(species[i]);

        I += pow(Zi,2) * Ci;
   }

   I.dimensions().reset(dimless);

   return 0.5 * I;

}

Foam::volScalarField
Foam::nonIdealST::gamma(const word& specie) const
{
    const scalar& Zi = speciesComposition().Z(specie);

    const volScalarField& I = this->I();
    
    volScalarField Gamma
    (
        IOobject
        (
            "Gamma", 
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero",dimless, gamma_)
    );   

    if(gammaVariation_ == true)
    {
        Gamma == ADH_*pow(Zi,2)*(pow(I, 0.5)/(1 + pow(I, 0.5))) - 0.3 * I; 
        Gamma.max(0);
        Gamma = exp(-Gamma);
    }

    return Gamma;
}

/****************************************************/
/*               Saturation Indexes                 */
/****************************************************/

Foam::volScalarField Foam::nonIdealST::IAP() const
{
    return volScalarField::null();
}


// Supersaturation ratio
Foam::volScalarField Foam::nonIdealST::SR() const
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
        volScalarField Gammai(this->gamma(species[i])); 

        tSR *=  pow((Gammai * Ci), rExp);
    }

    tSR /= Ksp_;

    tSR.max(0);

    return tSR;
}

// Supersaturation index
Foam::volScalarField Foam::nonIdealST::SI() const
{
    volScalarField tSI(this->SR());
    
    tSI.max(1.0);
    tSI = log10(tSI);

    return tSI;
}

// Absolute supersaturation 
Foam::volScalarField Foam::nonIdealST::S() const
{
    volScalarField tS(this->SR());

    tS = pow(Ksp_.value(), (1/nReactants_)) * (pow(tS, (1/nReactants_)) - 1);
    tS.max(0);
   
    // tS.dimensions().reset(dimensionSet(0,-3,0,0,1,0,0));

    return tS;
}

// Relative supersaturation
Foam::volScalarField Foam::nonIdealST::relS() const
{
    volScalarField tRelS(this->SR());

    tRelS = pow(tRelS, (1/nReactants_)) - 1;
    tRelS.max(0);

    return tRelS;
}

// Activity supersaturation or sqrt(Supersaturation ratio)
Foam::volScalarField Foam::nonIdealST::Sa() const
{
    volScalarField tSa(this->SR());

    tSa = gamma_ * pow(tSa, (1/nReactants_));

    return tSa;
}

void nonIdealST::correct()
{}


bool nonIdealST::read()
{
    return true;
}
void nonIdealST::printThermoParams()
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

        Info<<"Chemical Conversion = "<< max(C("MAP")).value()/(0.04 - pow(Ksp_.value(), 0.33)) << endl;


        Info<<"/*************************************************************/" << endl;
    }
}

void nonIdealST::updateConcentrations()
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
