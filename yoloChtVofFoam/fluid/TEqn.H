{
    

    const dimensionedScalar Cp1 = thermo.Cp1();
    const dimensionedScalar Cp2 = thermo.Cp2();

   
    rhoCp = alpha1*rho1*Cp1 + alpha2*rho2*Cp2;



    kappaEff = thermo.kappa() + rhoCp*turbulence.nut()/Prt;


    const surfaceScalarField rhoCpPhi
    (
        "rhoCpPhi",
        alphaPhi10*(rho1*Cp1 - rho2*Cp2) + phi*rho2*Cp2
    );

    fvScalarMatrix TEqn
    (
         fvm::ddt(rhoCp, T)
       + fvm::div(rhoCpPhi, T)
       - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
       - fvm::laplacian(kappaEff, T)
       + mixture.TSource()
    );


    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
