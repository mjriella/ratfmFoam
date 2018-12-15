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

#include "dispersedkEpsilon.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseTurbulenceModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dispersedkEpsilon, 0);
addToRunTimeSelectionTable(particleTurbulence,dispersedkEpsilon,dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dispersedkEpsilon::dispersedkEpsilon
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dragModel& drag1
)
:
    particleTurbulence(phase1,phase2,alpha1,alpha2,drag1),
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
    k_
    (
        IOobject
        (
            "k1",
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
            "epsilon1",
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
            "nut1",
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dispersedkEpsilon::~dispersedkEpsilon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool dispersedkEpsilon::read()
{
    // This could be better
    Cmu_.readIfPresent(coeffsDict());
    C1_.readIfPresent(coeffsDict());
    C2_.readIfPresent(coeffsDict());
    sigmaEps_.readIfPresent(coeffsDict());

    return true;
}


void dispersedkEpsilon::solve()
{
	particleTurbulence::solve();

    volVectorField U = phase1_.U();
    volVectorField Uopp = phase2_.U();
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(2*nut_*(tgradU() && dev(symm(tgradU()))));
    tgradU.clear();

	// kinetic viscosity
    volScalarField nuKT  =  
                phase1_.U().mesh().lookupObject<volScalarField>("nuKT");

    // drag term
    volScalarField beta = 2*(1-alpha1_)*drag1_.K(mag(U-Uopp))/rho1_;
   
    // coupling
    volScalarField k2 =  
                phase2_.U().mesh().lookupObject<volScalarField>("k2");
    volScalarField epsilon2  =  
                phase2_.U().mesh().lookupObject<volScalarField>("epsilon2");
    volScalarField kCoup = sqrt(k_*k2);
    volScalarField epsilonCoup = sqrt(epsilon_*epsilon2);
 
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi1_, epsilon_)
      - fvm::Sp(fvc::div(phi1_), epsilon_)
      - fvm::laplacian(nut_/sigmaEps_ + nuKT, epsilon_)
    ==
        C1_*G*epsilon_/k_ 
      - fvm::Sp(C2_*epsilon_/k_, epsilon_) 
      + beta*epsilonCoup
      - fvm::SuSp(beta, epsilon_)
    );
    epsEqn.relax();
    epsEqn.solve();
    bound(epsilon_, epsilonMin_);
    
    fvScalarMatrix kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi1_, k_)
      - fvm::Sp(fvc::div(phi1_), k_)
      - fvm::laplacian(nut_ + nuKT, k_)
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
    nut_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multiphaseTurbulenceModels
} // End namespace Foam

// ************************************************************************* //
