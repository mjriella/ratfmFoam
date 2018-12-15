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

#include "v2f.H"
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

defineTypeNameAndDebug(v2f, 0);
addToRunTimeSelectionTable(fluidTurbulence,v2f,dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<volScalarField> v2f::davidsonCorrectNut
(
    const tmp<volScalarField>& value
) const
{
    return min(CmuKEps_*sqr(k_)/epsilon_, value);
}


tmp<volScalarField> v2f::Ts() const
{
    return max(k_/epsilon_, 6.0*sqrt(nu2_/epsilon_));
}


tmp<volScalarField> v2f::Ls() const
{
    return CL_*max(pow(k_, 1.5)/epsilon_, Ceta_*pow025(pow3(nu2_)/epsilon_));
}

void v2f::correctEpsilon()
{
	tmp<volScalarField> nu_temp
	(
		new volScalarField
		(
    		IOobject
    		(
        		"nu_temp",
            	phase1_.U().time().timeName(),
            	phase1_.U().mesh(),
        		IOobject::NO_READ,
        		IOobject::NO_WRITE
    		),
    		phase1_.U().mesh(),
    		dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), nu2_.value())
		)
	);

	volScalarField& nu2 = nu_temp();

    const fvPatchList& patches = phase1_.U().mesh().boundary();
    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];
		const scalarField& nu2f = nu2.boundaryField()[patchi];

        if (isA<wallFvPatch>(p))
        {
			forAll(p, facei)
			{
				label faceCelli = p.faceCells()[facei];
				
				epsilon_[faceCelli] =
					2*k_[faceCelli]*nu2f[facei]/
							sqr(y_[patchi][facei]);
			}
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2f::v2f
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
            0.22
        )
    ),
    CmuKEps_
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
            1.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffsDict_,
            0.3
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffsDict_,
            0.23
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceta",
            coeffsDict_,
            70.0
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffsDict_,
            1.9
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            coeffsDict_,
            1.0
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
    ),
    v2_
    (
        IOobject
        (
            "v2",
            phase1.U().time().timeName(),
            phase1.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase1.U().mesh()
    ),
    f_
    (
        IOobject
        (
            "f2",
            phase1.U().time().timeName(),
            phase1.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase1.U().mesh()
    ),
    v2Min_(dimensionedScalar("v2Min", v2_.dimensions(), SMALL)),
    fMin_(dimensionedScalar("fMin", f_.dimensions(), 0.0))
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(v2_, v2Min_);
    bound(f_, fMin_);
    nut_ = davidsonCorrectNut(Cmu_*v2_*Ts());
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

v2f::~v2f()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool v2f::read()
{
    // This could be better
    Cmu_.readIfPresent(coeffsDict());
    CmuKEps_.readIfPresent(coeffsDict());
    C1_.readIfPresent(coeffsDict());
    C2_.readIfPresent(coeffsDict());
    CL_.readIfPresent(coeffsDict());
    Ceta_.readIfPresent(coeffsDict());
    Ceps2_.readIfPresent(coeffsDict());
    sigmaK_.readIfPresent(coeffsDict());
    sigmaEps_.readIfPresent(coeffsDict());

    return true;
}


void v2f::solve()
{
	fluidTurbulence::solve();

    // use N=6 so that f=0 at walls
    const dimensionedScalar N("N", dimless, 6.0);
    const volScalarField T(Ts());
    const volScalarField L2(type() + ".L2", sqr(Ls()));
    
    volScalarField alpha
    (
        "v2f::alpha",
        1.0/T*((C1_ - N)*v2_ - 2.0/3.0*k_*(C1_ - 1.0))
    );

    tmp<volScalarField> Ceps1 =
        1.4*(1.0 + 0.05*min(sqrt(k_/v2_), scalar(100.0)));

    volVectorField U = phase1_.U();
    volVectorField Uopp = phase2_.U();

    // Drag
    volScalarField beta = 2*alpha1_*drag1_.K(mag(Uopp-U))/rho2_;

    // Coupling
    volScalarField k1 =  
                phase2_.U().mesh().lookupObject<volScalarField>("k1");
    volScalarField epsilon1  =  
                phase2_.U().mesh().lookupObject<volScalarField>("epsilon1");
    volScalarField v1  =  
                phase2_.U().mesh().lookupObject<volScalarField>("v1");
	
    volScalarField kCoup = sqrt(k_*k1);
    volScalarField epsilonCoup = sqrt(epsilon_*epsilon1);
    volScalarField vCoup = sqrt(v2_*v1);

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(2*nut_*(tgradU() && dev(symm(tgradU()))));
    tgradU.clear();

	correctEpsilon(); // Enforce BC
    
    // turbulence kinetic energy dissipation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi1_, epsilon_)
      - fvm::Sp(fvc::div(phi1_), epsilon_)
      - fvm::laplacian(nut_/sigmaEps_ + nu2_, epsilon_)
    ==
        1.59*G/T // Ceps1_
      - fvm::Sp(Ceps2_/T, epsilon_) 
      + beta*epsilonCoup
      - fvm::SuSp(beta, epsilon_)
    );
    epsEqn().relax();
    epsEqn().solve();
    bound(epsilon_, epsilonMin_);
    
    // turbulence kinetic energy
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi1_, k_)
      - fvm::Sp(fvc::div(phi1_), k_)
      - fvm::laplacian(nut_/sigmaK_ + nu2_, k_)
    ==
        G
      - fvm::Sp(epsilon_/k_, k_)
      + beta*kCoup
      - fvm::SuSp(beta, k_)
    );
    kEqn().relax();
    kEqn().solve();
    bound(k_, kMin_);

    // Relaxation function equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2, f_)
      - 1.0/L2/k_*(alpha - C2_*G)
    );
    fEqn().relax();
    fEqn().solve();
    bound(f_, fMin_);

    // Normal component 
    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi1_, v2_)
      - fvm::laplacian(nut_/sigmaK_ + nu2_, v2_)
      ==
        min(k_*f_, C2_*G - alpha)
      - fvm::Sp(N*epsilon_/k_, v2_)
	  + beta*vCoup
	  - fvm::SuSp(beta, v2_)
    );
    v2Eqn().relax();
    v2Eqn().solve();
    bound(v2_, v2Min_);

    nut_ = davidsonCorrectNut(Cmu_*v2_*T);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multiphaseTurbulenceModels
} // End namespace Foam

// ************************************************************************* //
