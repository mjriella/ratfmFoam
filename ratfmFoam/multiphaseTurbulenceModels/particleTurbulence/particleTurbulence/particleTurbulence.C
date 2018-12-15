/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "particleTurbulence.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseTurbulenceModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(particleTurbulence, 0);
defineRunTimeSelectionTable(particleTurbulence, dictionary);
addToRunTimeSelectionTable(multiphaseTurbulenceModel, particleTurbulence, multiphaseTurbulenceModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleTurbulence::particleTurbulence
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dragModel& drag1
)
:
    multiphaseTurbulenceModel(phase1, phase2, alpha1, alpha2, drag1),
    multiphaseTurbulenceProperties_
    (
        IOobject
        (
            "multiphaseTurbulenceProperties",
            phase1_.U().time().constant(),
            phase1_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    coeffsDict_(multiphaseTurbulenceProperties_.subOrEmptyDict("coeffs")),
    kMin_("kMin", sqr(dimVelocity), SMALL),
    epsilonMin_("epsilonMin", kMin_.dimensions()/dimTime, SMALL)
{
    kMin_.readIfPresent(multiphaseTurbulenceProperties_);
    epsilonMin_.readIfPresent(multiphaseTurbulenceProperties_);
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<particleTurbulence> particleTurbulence::New
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
        ).lookup("particleTurbulence")
    );

    Info<< "Selecting particle turbulence model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RASModel::New"
            "("
                "const volVectorField&, "
                "const surfaceScalarField&, "
                "transportModel&, "
                "const word&"
            ")"
        )   << "Unknown RASModel type "
            << modelType << nl << nl
            << "Valid RASModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<particleTurbulence>
    (
        cstrIter()(phase1, phase2, alpha1, alpha2, drag1)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleTurbulence::~particleTurbulence()
{}


// * * * * * * * * * * * * * * * * Member Function   * * * * * * * * * * * * //
void particleTurbulence::solve()
{
    multiphaseTurbulenceModel::solve();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End multiphaseTurbulenceModels
} // End Foam
// ************************************************************************* //
