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

#include "MAPST.H"
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

defineTypeNameAndDebug(MAPST, 0);
addToRunTimeSelectionTable(speciesThermoModel, MAPST, speciesThermoModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MAPST::MAPST
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
    gamma_(coeffsDict_.lookupOrDefault<scalar>("gamma", 1)),
    ADH_(coeffsDict_.lookupOrDefault<scalar>("ADH", 1)),
   // kv_(readScalar(coeffsDict_.lookup("kv"))),
   // ka_(readScalar(coeffsDict_.lookup("ka"))),

    // Equilibrium constants
    equilibriumCoeffs_(coeffsDict_.subDict("equilibriumConstants")),
    kMgOH_
    (
        "kMgOH", 
        dimensionSet(0,-3,0,0,1,0,0), 
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKMgOH")))
    ),
    kMgHPO4_
    (
        "kMgHPO4", 
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKMgHPO4")))
    ),
    kMgH2PO4_
    (
        "kMgH2PO4",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKMgH2PO4")))
    ),
    kMgPO4_
    (
        "kMgPO4",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKMgPO4")))
    ),
    kH3PO4_
    (
        "kH3PO4",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKH3PO4")))
    ),
    kH2PO4_
    (
        "kH2PO4",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKH2PO4")))
    ),
    kHPO4_
    (
        "kHPO4",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKHPO4")))
    ),
    kNH4_
    (
        "kNH4",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKNH4")))
    ),
    kH2O_
    (
        "kH2O",
        //dimless, 
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,readScalar(equilibriumCoeffs_.lookup("pKH2O")))
    ),

    Ksp_
    (
        "Ksp",
        dimensionSet(0,-9,0,0,3,0,0), 
        pow(10, readScalar(equilibriumCoeffs_.lookup("pKsp")))
    ),

    // Ion concentration fields

    C_(this->speciesComposition().species().size()),
    Ct_(this->speciesComposition().species().size())
        
{

    // Populate C_ and Ct_
   
    Info<<"DEEP CODE"<<endl;
    
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

    // Print equilibrium constants
    Info<<"kHPO4 = "  << kHPO4_    << "\n"
        <<"kH2PO4 = " << kH2PO4_   << "\n"
        <<"kH3PO4 = " << kH3PO4_   << "\n"
        <<"kMGH2PO4= "<< kMgH2PO4_ << "\n"
        <<"kMGHPO4 = "<< kMgHPO4_  << "\n"
        <<"kMGPO4= "  << kMgPO4_   << "\n"
        <<"kMGOF= "   << kMgOH_    << "\n"
        <<"kNH4 = "   << kNH4_     << "\n"
        <<"Ksp = "    << Ksp_      << "\n"
        <<"kH2O= "    << kH2O_     << endl;

    Info<<"ADH = "<< ADH_ << endl;

}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<MAPST> MAPST::New
(
    volVectorField& U,
    volScalarField& rho,
    dictionary& speciesDict
)
{
    return autoPtr<MAPST>
    (
        new MAPST(U, rho, speciesDict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::MAPST::R
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


/****************************************************/
/*             Speciation constants                 */
/****************************************************/

Foam::dimensionedScalar Foam::MAPST::Hconc() const 
{
    dimensionedScalar Hc_
    (
        "Hc_",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10,-pH_)   
    );
        
    return Hc_;
}

Foam::dimensionedScalar Foam::MAPST::OHconc() const
{
    dimensionedScalar OHc_
    (   
        "OHc_",
        dimensionSet(0,-3,0,0,1,0,0),
        pow(10, -(14-pH_))
    );

    return OHc_;
}       

Foam::volScalarField Foam::MAPST::K0() const
{
    // H+ concentration
    const dimensionedScalar& Hc_ = this->Hconc();

    volScalarField tK0
    (
        IOobject
        (
            "K0",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        (
            (pow(Hc_, 3)/(kH3PO4_*kH2PO4_*kHPO4_)) +
            (pow(Hc_, 2)/(kH2PO4_*kHPO4_)) +
            (Hc_/kHPO4_) +
            1
        )
    );
    return tK0;
}

Foam::volScalarField Foam::MAPST::K1() const
{
    // H+ concentration
    const dimensionedScalar& Hc_ = this->Hconc();

    volScalarField tK1
    (
        IOobject
        (
            "K1",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        (
            (pow(Hc_, 2)/(kMgH2PO4_*kH2PO4_*kHPO4_)) +
            (Hc_/(kMgHPO4_*kHPO4_)) +
            1/kMgPO4_
        )
    );
    return tK1;
}

Foam::volScalarField Foam::MAPST::K2() const
{
    // OH- concentration
    const dimensionedScalar& OHc_ = this->OHconc();

    volScalarField tK2
    (
       IOobject
       (
           "K2",
           this->mesh().time().timeName(),
           this->mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE,
           false
       ),
       this->mesh(),
       (OHc_/kMgOH_) +1
    );
    
    return tK2;
}

Foam::volScalarField Foam::MAPST::K3() const
{
    // H+ concentration
    const dimensionedScalar& Hc_ = this->Hconc();

    volScalarField tK3
    (
        IOobject
        (
            "K3",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        (kNH4_/Hc_) + 1
    );

    return tK3;
}

/****************************************************/
/*               Ion concentrations                 */
/****************************************************/

Foam::volScalarField Foam::MAPST::C_MG() const
{
    // Refs to Ks
    const volScalarField& K0 = this->K0();
    const volScalarField& K1 = this->K1();
    const volScalarField& K2 = this->K2();

    // Mass fractions (C_Ts in Ye et al 2018)
    const volScalarField& Ct_PO4 = this->Ct("PO4");
    const volScalarField& Ct_MG = this->Ct("MG");
    
    return (
                (
                  K1*Ct_MG - K0*K2 - K1*Ct_PO4
                + sqrt( 
                        pow((K1*Ct_MG - K0*K2 - K1*Ct_PO4), 2) 
                      + 4*K0*K1*K2*Ct_MG)
                )/ (2*K1*K2)
           );

}


Foam::volScalarField Foam::MAPST::C_NH4() const
{
    // Total concentration of NH4
    const volScalarField& Ct_NH4 = this->Ct("NH4");

    return Ct_NH4/this->K3();
}

Foam::volScalarField Foam::MAPST::C_PO4() const
{
    const volScalarField& K0 = this->K0();
    const volScalarField& K1 = this->K1();

    // Total concentration
    const volScalarField& Ct_PO4 = this->Ct("PO4");

    return Ct_PO4/(K0+(K1*this->C_MG()));
}

Foam::volScalarField Foam::MAPST::C_MAP() const
{
    return this->Ct("MAP");
}

/****************************************************/
/*             Ionic activity functions             */
/****************************************************/

Foam::volScalarField
Foam::MAPST::I() const
{

   volScalarField I
   (
        IOobject
        (
            "I", 
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero",dimMoles/dimVolume,0)
   );

   const speciesTable& species = speciesComposition().species(); 

   forAll(species, i)
   {
        const scalar& Zi = this->speciesComposition().Z(i);
        const volScalarField& Ci = this->Ct(species[i]);

        I += pow(Zi,2) * Ci;
   }

  // Info<<"min(I) = "<< min(I).value()
  //     <<" max(I) = " << max(I).value() << endl;

   I.dimensions().reset(dimless);

   return 0.5 * I;

}

Foam::volScalarField
Foam::MAPST::gamma(const word& specie) const
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

Foam::volScalarField Foam::MAPST::IAP() const
{
    const volScalarField& C_MG = this->C("MG");
    const volScalarField& C_PO4 = this->C("PO4");
    const volScalarField& C_NH4 = this->C("NH4");

 //   Info<<"max(Gamma(MG)) = "<< max(this->gamma("MG")).value() << "\n"
 //       <<"max(Gamma(PO4)) = "<< max(this->gamma("PO4")).value() << "\n"
 //       <<"max(Gamma(NH4)) = "<< max(this->gamma("NH4")).value() << endl;
   
    volScalarField tIAP(  this->gamma("MG")  * C_MG 
                        * this->gamma("PO4") * C_PO4 
                        * this->gamma("NH4") * C_NH4 );
    
    return tIAP;

}

// Supersaturation ratio
Foam::volScalarField Foam::MAPST::SR() const
{
    volScalarField tSR(this->IAP()/Ksp_);
    tSR.max(0);

    return tSR;
}

// Supersaturation index
Foam::volScalarField Foam::MAPST::SI() const
{
    volScalarField tSI(this->SR());
    tSI.max(1.0);
    tSI = log10(tSI);

    return tSI;
}

// Absolute supersaturation 
Foam::volScalarField Foam::MAPST::S() const
{
    volScalarField tS(pow(this->IAP(), 0.33) - pow(Ksp_, 0.33));
    tS.max(0);
   
    tS.dimensions().reset(dimless);

    return tS;
}

// Relative supersaturation
Foam::volScalarField Foam::MAPST::relS() const
{
    volScalarField tRelS(this->SR());

    tRelS = pow(tRelS, 0.33) - 1;
    tRelS.max(0);

    return tRelS;
}

// Activity supersaturation or cbrt(Supersaturation ratio)
Foam::volScalarField Foam::MAPST::Sa() const
{
    volScalarField tSa(this->SR());

    tSa = pow(tSa, 0.33);

    return tSa;
}

/***********************************************************/

void MAPST::correct()
{}


bool MAPST::read()
{
    return true;
}

void Foam::MAPST::printThermoParams()
{

    if(printParams_)
    {
        
        Info<<"\n/*************************************************************/ \n" 
            <<"/**********   Supersaturation parameters output: ***************/ \n"
            <<"/***************************************************************/ \n"

            <<"Saturation ratio min(SR) = "<< min(this->SR()).value()
            <<" max(SR) = "<< max(this->SR()).value() << "\n"

            <<"Saturation Index min(SI) = "<< min(this->SI()).value()
            <<" max(SI) = "<< max(this->SI()).value() << "\n"

            <<"Absolute saturation min(S) = "<< min(this->S()).value()
            <<" max(S) = "<< max(this->S()).value() << "\n"

            <<"Activity saturation ratio min(Sa) = "<< min(this->Sa()).value()
            <<" max(Sa) = "<< max(this->Sa()).value() <<  "\n"

            <<"Relative supersaturation min(relS) = "<< min(this->relS()).value()
            <<" max(relS) = "<< max(this->relS()).value() << "\n\n" 
            << "***********  Concentrations  ************\n" << endl;

        scalar reqH2O = 0;

        forAll(this->C(),i)
        {
            volScalarField& Ci = this->C(i);
            volScalarField& Cti = this->Ct(i);
            dimensionedScalar& Wi = this->speciesComposition().W(i);

            const word& specie = this->speciesComposition().Y(i).name();
            const scalar reqC = this->speciesDict().subDict(specie).lookupOrDefault<scalar>("requiredC", 0.0);

            if(specie != "H2O")
            {
                Info<<"Ion Concentration: "
                    <<"min({"       << specie << "}) = "  << min(Ci).value() 
                    <<" max({"      << specie << "}) = "  << max(Ci).value() << "\n"
                    << "Total Concentration: "
                    <<"min(["      << specie <<"]) = "  << min(Cti).value() 
                    <<" max(["     << specie << "]) = " << max(Cti).value() << "\n" 
                    <<"Required(Y_"<< specie << ") = "  << reqC*Wi.value()/max(this->rho()).value() << endl; 

                reqH2O += reqC*Wi.value()/max(this->rho()).value();

            } else {
                Info<<"min("<< specie << ") = "<< min(Ci).value() 
                    <<" max("<< specie << ") = "<< max(Ci).value() << "\n"
                    <<"Required(Y_"<< specie << ") = " << 1 - reqH2O << endl;  
            }
        }
        
        Info<<"Chemical Conversion = "<< max(Ct("MAP")).value()/(0.04 - pow(Ksp_.value(), 0.33)) << endl;

        
        Info <<"/*************************************************************/\n" << endl;
    }
}

void Foam::MAPST::updateConcentrations()
{

    forAll(this->Ct(), i)
    {  
        // Update total concentration fields from mass fractions
        volScalarField& Cti = this->Ct(i);
        volScalarField& Yi = this->speciesComposition().Y(i);
        dimensionedScalar& Wi = this->speciesComposition().W(i);
        volScalarField& rho = this->rho();

        Cti = rho * Yi / Wi;
        
        // Update ion concentration fields (might be redundant)
        // TODO SPAGHETTI CODE have to do that because naming
        //      and because of equations for ion concentrations
        //      fix that when clean up is being done
        volScalarField& Ci = this->C(i);
        
        if(Ci.name() == "C_PO4")
        {
            Ci = this->C_PO4();
        } else if (Ci.name() == "C_NH4")
        {
            Ci = this->C_NH4();
        } else if (Ci.name() == "C_MG")
        {
            Ci = this->C_MG();
        } else if (Ci.name() == "C_MAP")
        {
            Ci = this->C_MAP();
        } else {
            Ci = rho * Yi / Wi;  // this is for water ion field TODO Remove that field
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
