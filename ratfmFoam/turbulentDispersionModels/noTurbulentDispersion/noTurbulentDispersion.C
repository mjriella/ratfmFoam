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

#include "noTurbulentDispersion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(noTurbulentDispersion, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        noTurbulentDispersion,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::noTurbulentDispersion::noTurbulentDispersion
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

Foam::turbulentDispersionModels::noTurbulentDispersion::
~noTurbulentDispersion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::noTurbulentDispersion::D() const
{
    
	return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "zero",
                    alpha1_.mesh().time().timeName(),
                    alpha1_.mesh()
                ),
                alpha1_.mesh(),
                dimensionedScalar
                (
                    "zero",
                    dimD,
                    0
                )
            )
        );
}
/*
Foam::tmp<Foam::volVectorField>
Foam::turbulentDispersionModels::noTurbulentDispersion::D() const
{
    const fvMesh& mesh = phase1_.mesh();
    return
        tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "zero",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, 1, -1, 0, 0),
                    vector::zero
                )
            )
        );
}
*/

// ************************************************************************* //
