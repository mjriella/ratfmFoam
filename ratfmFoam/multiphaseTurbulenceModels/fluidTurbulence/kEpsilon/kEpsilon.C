/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "kEpsilon.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseTurbulenceModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilon, 0);
addToRunTimeSelectionTable(fluidTurbulence,kEpsilon,dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilon::kEpsilon
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dragModel& drag1
)
:
    fluidTurbulence(phase1,phase2,alpha1,alpha2,drag1),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffsDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffsDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffsDict_,
            1.92
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffsDict_,
            1.3
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffsDict_,
            0.41
        )
    ),
    E_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffsDict_,
            9.8
        )
    ),
    k_
    (
        IOobject
        (
            "k2",
            phase1.U().time().timeName(),
            phase1.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase1.U().mesh()
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon2",
            phase1.U().time().timeName(),
            phase1.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase1.U().mesh()
    ),
    nut_
    (
        IOobject
        (
            "nut2",
            phase1.U().time().timeName(),
            phase1.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase1.U().mesh()
    )
{
    bound(epsilon_, epsilonMin_);
    bound(k_, kMin_);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void kEpsilon::correctEpsilon()
{
    labelList cellBoundaryFaceCount(epsilon_.size(), 0);

    scalar Cmu25 = ::pow(Cmu_.value(), 0.25);
    scalar Cmu75 = ::pow(Cmu_.value(), 0.75);
    scalar kappa = kappa_.value();
    scalar nu2 = nu2_.value();
	
	// TODO this should be changed as its a bit quick and dirty
    volVectorField U = phase1_.U();
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(2*nut_*(tgradU() && dev(symm(tgradU()))));
    tgradU.clear();

    const fvPatchList& patches = phase1_.U().mesh().boundary();

    //- Initialise the near-wall P field to zero
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                epsilon_[faceCelli] = 0.0;
                G[faceCelli] = 0.0;
            }
        }
    }

    // Enhanced wall treatment 
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const scalarField& nut2w = nut_.boundaryField()[patchi];
            scalarField magFaceGradU(mag(phase1_.U().boundaryField()[patchi].snGrad()));

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];
                scalar yPlus = Cmu25*y_[patchi][facei]*::sqrt(k_[faceCelli])/nu2;
                scalar yPlusLam = 11.0;
                cellBoundaryFaceCount[faceCelli]++;

                if (yPlus > yPlusLam)
                {
                	epsilon_[faceCelli] +=
                     	Cmu75*::pow(k_[faceCelli], 1.5)
                    	/(kappa*y_[patchi][facei]);

                	G[faceCelli] +=
                    	(nut2w[facei] + nu2)*magFaceGradU[facei]
                   		*Cmu25*::sqrt(k_[faceCelli])
                   		/(kappa*y_[patchi][facei]);
                }
                else
                {
                    epsilon_[faceCelli] += 2.0*k_[faceCelli]*nu2/
                                    sqr(y_[patchi][facei]);
                }
            }
        }
    }
}

void kEpsilon::correctNut()
{
    scalar Cmu25 = ::pow(Cmu_.value(), 0.25);
    scalar kappa = kappa_.value();
    scalar E = E_.value();
    scalar nu2 = nu2_.value();

    const fvPatchList& patches = phase1_.U().mesh().boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nut_.boundaryField()[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                // calculate yPlus
                scalar yPlus =
                    Cmu25*y_[patchi][facei]*::sqrt(k_[faceCelli])/nu2;

                if (yPlus > 11.6)
                {
                    nutw[facei] = nu2*(yPlus*kappa/::log(E*yPlus) -1);
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kEpsilon::~kEpsilon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool kEpsilon::read()
{
    // This could be better
    Cmu_.readIfPresent(coeffsDict());
    C1_.readIfPresent(coeffsDict());
    C2_.readIfPresent(coeffsDict());
    sigmaEps_.readIfPresent(coeffsDict());

    return true;
}


void kEpsilon::solve()
{
	fluidTurbulence::solve();

    volVectorField U = phase1_.U();
    volVectorField Uopp = phase2_.U();
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(2*nut_*(tgradU() && dev(symm(tgradU()))));
    tgradU.clear();

    // drag term
    volScalarField beta = 2*alpha1_*drag1_.K(mag(Uopp-U))/rho2_;
   
    // coupling
    volScalarField k1 =  
                phase2_.U().mesh().lookupObject<volScalarField>("k1");
    volScalarField epsilon1  =  
                phase2_.U().mesh().lookupObject<volScalarField>("epsilon1");
    volScalarField kCoup = sqrt(k_*k1);
    volScalarField epsilonCoup = sqrt(epsilon_*epsilon1);

	correctEpsilon(); // Wall func for epsilon
    
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi1_, epsilon_)
      - fvm::Sp(fvc::div(phi1_), epsilon_)
      - fvm::laplacian(nut_/sigmaEps_ + nu2_, epsilon_)
    ==
        C1_*G*epsilon_/k_ 
      - fvm::Sp(C2_*epsilon_/k_, epsilon_) 
      + beta*epsilonCoup
      - fvm::SuSp(beta, epsilon_)
    );
    epsEqn.relax();
	epsEqn.boundaryManipulate(epsilon_.boundaryField());
    epsEqn.solve();
    bound(epsilon_, epsilonMin_);
    
    fvScalarMatrix kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi1_, k_)
      - fvm::Sp(fvc::div(phi1_), k_)
      - fvm::laplacian(nut_ + nu2_, k_)
    ==
        G
      - fvm::Sp(epsilon_/k_, k_)
      + beta*kCoup
      - fvm::SuSp(beta, k_)
    );
    kEqn.relax();
    kEqn.solve();
    bound(k_, kMin_);

    nut_ = Cmu_*sqr(k_)/epsilon_;
	correctNut();
    nut_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multiphaseTurbulenceModels
} // End namespace Foam

// ************************************************************************* //
