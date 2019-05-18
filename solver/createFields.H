Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading incompressibleTwoPhaseInteractingMixture\n" << endl;
incompressibleTwoPhaseInteractingMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());

// Mixture density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.rho()
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField rhod
(
    IOobject
    (
        "rhod",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mixture.rhod()
);

// Relative Velocity
autoPtr<relativeVelocityModel> UdmModelPtr
(
    relativeVelocityModel::New
    (
        mixture,
        mixture
    )
);


// Construct compressible turbulence model
autoPtr
<
    CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
> turbulence
(
    CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
    ::New(rho, U, rhoPhi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

#include "createMRF.H"

// Species Part

    // Instantiate precipitation model/class
    autoPtr<Foam::precipitationModel> precipitationModel
    (
        Foam::precipitationModel::New(U, rho, phi)
    );
    
  //  const word& precipitationModelName = precipitationModel->precipitationModelName();
  //  const word& alphaPrecipitationModel = precipitationModel->alphaPrecipitationModelName();

    // Species reference
    speciesThermoModel& speciesThermo= precipitationModel->speciesThermo();
    speciesMixture& composition = speciesThermo.speciesComposition();

    
    // Species fields
    PtrList<volScalarField>& Y = composition.Y();
    
    word inertSpecie(precipitationModel->speciesDict().lookup("inertSpecie"));

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }


    // Thermodynamic saturation indicators
    speciesThermo.updateConcentrations();

    volScalarField SI
    (
        IOobject
        (
            "SI",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),        
        speciesThermo.SI()
    );    

    volScalarField Sa
    (
        IOobject
        (
            "Sa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),        
        speciesThermo.Sa()
    ); 

    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),        
        speciesThermo.S()
    );

    
    // Saueter diameter field
    volScalarField& d43 = mixture.d43();
    volScalarField& d32 = mixture.d43();

    volVectorField Vt
    (
        IOobject
        (
            "Vt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        g*d43*d43*(mixture.rhod() - mixture.rhoc())/(18 * mixture.mu())
    );

    speciesThermo.printThermoParams();
   
    const IOdictionary& pbeDict = precipitationModel->populationBalanceProperties();
    dimensionedScalar Cg = 
           pbeDict.subDict("univariateCoeffs").subDict("growthModel").lookup("Cg");

    Cg.dimensions().reset(dimensionSet(0,1,-1,0,0,0,0));    

    // ALPHA RATIO HOT FIX TODO 
    volScalarField alphaRatio
    (
        IOobject
        (
            "alphaRatio",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1/(alpha1 + 1e-15),
        zeroGradientFvPatchScalarField::typeName
    );



