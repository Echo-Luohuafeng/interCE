// Initialise fluid field pointer lists
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<twoPhaseMixtureEThermo> thermoFluid(fluidRegions.size());
PtrList<temperaturePhaseChangeTwoPhaseMixture> mixtureFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
PtrList<Yolo::interfaceProperties> interfaceFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());
PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());
PtrList<volScalarField> kappaEffFluid(fluidRegions.size());
PtrList<volScalarField> rhoCpFluid(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size(), 0);
List<scalar> pRefValueFluid(fluidRegions.size(), 0.0);

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

List<bool> frozenFlowFluid(fluidRegions.size(), false);


forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid region properties: " << fluidRegions[i].name() << endl;

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject("p_rgh", runTime.timeName(), fluidRegions[i], IOobject::MUST_READ, IOobject::AUTO_WRITE),
            fluidRegions[i]
        )
    );

    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject("U", runTime.timeName(), fluidRegions[i], IOobject::MUST_READ, IOobject::AUTO_WRITE),
            fluidRegions[i]
        )
    );

    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject("phi", runTime.timeName(), fluidRegions[i], IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
            fvc::flux(UFluid[i])
        )
    );

    Info<< "    Creating twoPhaseMixtureEThermo" << endl;
    thermoFluid.set
    (
        i,
        new twoPhaseMixtureEThermo(UFluid[i], phiFluid[i])
    );

    Info<< "    Creating temperaturePhaseChangeTwoPhaseMixture" << endl;
    mixtureFluid.set
    (
        i,
        temperaturePhaseChangeTwoPhaseMixture::New(thermoFluid[i], fluidRegions[i]).ptr()
    );

    const volScalarField& alpha1 = thermoFluid[i].alpha1();
    const volScalarField& alpha2 = thermoFluid[i].alpha2();
    const dimensionedScalar& rho1 = thermoFluid[i].rho1();
    const dimensionedScalar& rho2 = thermoFluid[i].rho2();

    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject("rho", runTime.timeName(), fluidRegions[i], IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
            alpha1*rho1 + alpha2*rho2
        )
    );
    rhoFluid[i].oldTime();

    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject("rhoPhi", runTime.timeName(), fluidRegions[i], IOobject::NO_READ, IOobject::NO_WRITE),
            fvc::interpolate(rhoFluid[i])*phiFluid[i]
        )
    );

    interfaceFluid.set
    (
        i,
        new Yolo::interfaceProperties
        (
            alpha1,
            UFluid[i],
            thermoFluid[i].transportPropertiesDict()
        )
    );

    turbulenceFluid.set
    (
        i,
        incompressible::turbulenceModel::New(UFluid[i], phiFluid[i], thermoFluid[i]).ptr()
    );

    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject("hRef", runTime.constant(), fluidRegions[i], IOobject::READ_IF_PRESENT, IOobject::NO_WRITE),
            dimensionedScalar(word::null, dimLength, 0.0)
        )
    );

    dimensionedScalar ghRef = (mag(g.value()) > SMALL)
        ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
        : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0.0);

    ghFluid.set(i, new volScalarField("gh", (g & fluidRegions[i].C()) - ghRef));
    ghfFluid.set(i, new surfaceScalarField("ghf", (g & fluidRegions[i].Cf()) - ghRef));

    volScalarField& p = thermoFluid[i].p();
    setRefCell(p, p_rghFluid[i], fluidRegions[i].solutionDict().subDict("PIMPLE"), pRefCellFluid[i], pRefValueFluid[i]);

    if (p_rghFluid[i].needReference())
    {
        p += dimensionedScalar("p", p.dimensions(), pRefValueFluid[i] - getRefCellValue(p, pRefCellFluid[i]));
        p_rghFluid[i] = p - rhoFluid[i]*ghFluid[i];
    }

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
    fluidRegions[i].setFluxRequired(alpha1.name());

    MRFfluid.set(i, new IOMRFZoneList(fluidRegions[i]));
    fluidFvOptions.set(i, new fv::options(fluidRegions[i]));

    PrtFluid.set(i, new dimensionedScalar("Prt", dimless, thermoFluid[i].transportPropertiesDict()));

    kappaEffFluid.set
    (
        i,
        new volScalarField
        (
            IOobject("kappaEff", runTime.timeName(), fluidRegions[i], IOobject::NO_READ, IOobject::NO_WRITE),
            thermoFluid[i].kappa()
        )
    );

    rhoCpFluid.set
    (
        i,
        new volScalarField
        (
            IOobject("rhoCp", runTime.timeName(), fluidRegions[i], IOobject::NO_READ, IOobject::NO_WRITE),
            rhoFluid[i]*thermoFluid[i].Cp()
        )
    );
    rhoCpFluid[i].oldTime();
}
