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

#include "Zuber.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(Zuber, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        Zuber,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Zuber::Zuber
(
    const dictionary& interfaceDict,
    const volScalarField& alpha1,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    turbulentDispersionModel(interfaceDict, alpha1, phase1, phase2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Zuber::
~Zuber()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::Zuber::D() const
{
    const volScalarField& k2 = alpha1_.mesh().lookupObject<volScalarField> ("k2");
    const volScalarField& k1 = alpha1_.mesh().lookupObject<volScalarField> ("k1");
    const volScalarField& beta = alpha1_.mesh().lookupObject<volScalarField> ("beta");
    const volScalarField& nut2 = alpha1_.mesh().lookupObject<volScalarField> ("nut2");

    volScalarField alpha2 = 1.0 - alpha1_;
    volScalarField Scfp = sqrt(k2/k1);

    return  beta*nut2/*fvc::grad(alpha1_)*// // move grad alpha to p
            (
                max(alpha1_, scalar(1e-6))*
                max(alpha2, scalar(1e-6))*
                Scfp*phase1_.rho()
            );
}




// ************************************************************************* //
