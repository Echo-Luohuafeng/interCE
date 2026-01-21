twoPhaseMixtureEThermo& thermo = thermoFluid[i];

temperaturePhaseChangeTwoPhaseMixture& mixture = mixtureFluid[i];

volScalarField& alpha1 = thermo.alpha1();
volScalarField& alpha2 = thermo.alpha2();

const dimensionedScalar& rho1 = thermo.rho1();
const dimensionedScalar& rho2 = thermo.rho2();

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
surfaceScalarField& rhoPhi = rhoPhiFluid[i];

Yolo::interfaceProperties& interface = interfaceFluid[i];

incompressible::turbulenceModel& turbulence = turbulenceFluid[i];

volScalarField& p = thermo.p();
volScalarField& p_rgh = p_rghFluid[i];

const volScalarField& gh = ghFluid[i];
const surfaceScalarField& ghf = ghfFluid[i];

dimensionedScalar& Prt = PrtFluid[i];
volScalarField& kappaEff = kappaEffFluid[i];
volScalarField& rhoCp = rhoCpFluid[i];

IOMRFZoneList& MRF = MRFfluid[i];
fv::options& fvOptions = fluidFvOptions[i];

bool frozenFlow = false;
fluidRegions[i].solutionDict().subDict("PIMPLE").readIfPresent("frozenFlow", frozenFlow);

const label pRefCell = pRefCellFluid[i];
const scalar pRefValue = pRefValueFluid[i];

volScalarField& T = thermo.T();



IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

const bool alphaRestart =
    alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

if (alphaRestart)
{
    Info << "Restarting alpha" << endl;
}

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);
