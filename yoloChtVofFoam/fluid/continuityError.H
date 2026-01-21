
// volScalarField contErr(fvc::ddt(rho) + fvc::div(rhoPhi));


// Pair<tmp<volScalarField>> vDot = mixture->vDot();
// const volScalarField& vDotc = vDot[0]();
// const volScalarField& vDotv = vDot[1]();



volScalarField contErr(fvc::div(phi) - vDotc - vDotv );

scalar sumLocalContErr = runTime.deltaTValue()*
    mag(contErr)().weightedAverage(mesh.V()).value();

scalar globalContErr = runTime.deltaTValue()*
    contErr.weightedAverage(mesh.V()).value();

cumulativeContErr += globalContErr;

Info<< "time step continuity errors : sum local = " << sumLocalContErr
    << ", global = " << globalContErr
    << ", cumulative = " << cumulativeContErr
    << endl;
