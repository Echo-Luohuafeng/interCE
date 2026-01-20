/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "yoloInterfaceProperties.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::Yolo::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::Yolo::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    tmp<volVectorField> tgradAlpha;
    if (nAlphaSmoothCurvature_ < 1)
    {
        tgradAlpha = fvc::grad(alpha1_, "nHat");
    }
    else
    {
        // Smooth interface curvature to reduce spurious currents
        auto talpha1L = tmp<volScalarField>::New(alpha1_);
        auto& alpha1L = talpha1L.ref();

        for (int i = 0; i < nAlphaSmoothCurvature_; ++i)
        {
            alpha1L = fvc::average(fvc::interpolate(alpha1L));
        }

        tgradAlpha = fvc::grad(talpha1L, "nHat");
    }

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(tgradAlpha));

    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryFieldRef()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}

// SXD add : Compute the correction factor to surface tension force
Foam::tmp<Foam::volScalarField>
Foam::Yolo::interfaceProperties::cFact() const
{
    tmp<volScalarField> tCorr
    (
        new volScalarField
        (
            IOobject
            (
                "rhoCorr",
                alpha1_.mesh().time().timeName(),
                alpha1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            alpha1_.mesh(),
            dimensionedScalar("corr",dimless,1.)
        )
    );

    volScalarField& corr = tCorr.ref();

    if(rhoCorrection_)
    {
        wordList phases(transportPropertiesDict_.lookup("phases"));

        dimensionedScalar rhol
        (
            "rho",
            dimDensity,
            transportPropertiesDict_.subDict(phases[0]).get<scalar>("rho")
        );

        dimensionedScalar rhov
        (
            "rho",
            dimDensity,
            transportPropertiesDict_.subDict(phases[1]).get<scalar>("rho")
        );


        const  volScalarField& rho
        (
            alpha1_.mesh().lookupObject<volScalarField>("rho")
        );


        dimensionedScalar rhoMean((rhol+rhov)/2.0);


        corr = rho/rhoMean;
    }

    return tCorr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Yolo::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    nAlphaSmoothCurvature_
    (
        alpha1.mesh().solverDict(alpha1.name()).
            getOrDefault<int>("nAlphaSmoothCurvature", 0)
    ),
    cAlpha_
    (
        alpha1.mesh().solverDict(alpha1.name()).get<scalar>("cAlpha")
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(alpha1.mesh().V()))
    ),

    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, Zero)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),

    rhoCorrection_
    (
        alpha1.mesh().solverDict(alpha1.name()).
            getOrDefault<Switch>("rhoCorrection", false)
    )
{
    if(rhoCorrection_)
    {
        Info<<"interfaceProperties: Using density correction when computing"
            <<" the surface tension force."
            <<endl;
    }

    Info << "using alpha smooth : nAlphaSmoothCurvature is set " << nAlphaSmoothCurvature_<<endl;
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::Yolo::interfaceProperties::sigmaK() const
{
    return cFact()*(sigmaPtr_->sigma()*K_);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::Yolo::interfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::Yolo::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


void Foam::Yolo::interfaceProperties::correct()
{
    calculateK();
}

bool Foam::Yolo::interfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}




// ************************************************************************* //
