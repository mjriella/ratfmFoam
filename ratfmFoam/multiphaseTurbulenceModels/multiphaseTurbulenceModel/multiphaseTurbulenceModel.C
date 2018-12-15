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

#include "multiphaseTurbulenceModel.H"
#include "volFields.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(multiphaseTurbulenceModel, 0);
	defineRunTimeSelectionTable(multiphaseTurbulenceModel, multiphaseTurbulenceModel);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseTurbulenceModel::multiphaseTurbulenceModel
(
	const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dragModel& drag1
)
:
    phase1_(phase1),
	phase2_(phase2),
    U1_(phase1.U()),
    U2_(phase2.U()),
    alpha1_(alpha1),
	alpha2_(alpha2),
    phi1_(phase1.phi()),
    phi2_(phase2.phi()),
    drag1_(drag1),
    rho1_(phase1.rho()),
    rho2_(phase2.rho()),
    da_(phase1.d()),
    nu1_(phase1.nu()),
    nu2_(phase2.nu()),
    y_(phase1.U().mesh())
{}


// * * * * * * * * * * * * * * * * Selector    * * * * * * * * * * * * * * * //
Foam::autoPtr<Foam::multiphaseTurbulenceModel> 
Foam::multiphaseTurbulenceModel::New
(
	const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dragModel& drag1
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "multiphaseTurbulenceProperties",
                phase1.U().time().constant(),
                phase1.U().db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("multiphaseTurbulence" + phase1.name())
    );

    Info<< "Selecting multiphase turbulence model type " 
										<< modelType << endl;

    multiphaseTurbulenceModelConstructorTable::iterator cstrIter =
        multiphaseTurbulenceModelConstructorTablePtr_->find(modelType);

    if (cstrIter == multiphaseTurbulenceModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "multiphaseTurbulenceModel::New"
            "("
                "const phaseModel&, "
                "const phaseModel&, "
                "const volScalarField&, "
                "const volScalarField&,"
				"const dragModel&"
            ")"
        )   << "Unknown multiphaseTurbulenceModel type "
            << modelType << nl << nl
            << "Valid multiphaseTurbulenceModel types:" << endl
            << multiphaseTurbulenceModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<multiphaseTurbulenceModel>
    (
        cstrIter()(phase1, phase2, alpha1, alpha2, drag1)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseTurbulenceModel::~multiphaseTurbulenceModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseTurbulenceModel::solve()
{
    // polymorphism access through here
    // gets overloaded through derived class
}


// ************************************************************************* /}/
