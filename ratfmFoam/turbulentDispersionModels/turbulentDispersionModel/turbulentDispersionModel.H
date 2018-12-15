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

Class
    Foam::turbulentDispersionModel

Description

SourceFiles
    turbulentDispersionModel.C
    newTurbulentDispersionModel.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDispersionModel_H
#define turbulentDispersionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class turbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class turbulentDispersionModel
{
protected:

    // Protected data

        const dictionary& interfaceDict_;
        const volScalarField& alpha1_;
        const phaseModel& phase1_;
        const phaseModel& phase2_;


public:

    //- Runtime type information
    TypeName("turbulentDispersionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            turbulentDispersionModel,
            dictionary,
            (
                const dictionary& interfaceDict,
                const volScalarField& alpha1,
                const phaseModel& phase1,
                const phaseModel& phase2
            ),
            (interfaceDict, alpha1, phase1, phase2)
        );

    // Static data members
        
        //- Force the dimensions
        static const dimensionSet dimD;

    // Constructors

        turbulentDispersionModel
        (
            const dictionary& interfaceDict,
            const volScalarField& alpha1,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    //- Destructor
    virtual ~turbulentDispersionModel();


    // Selectors

        static autoPtr<turbulentDispersionModel> New
        (
            const dictionary& interfaceDict,
            const volScalarField& alpha1,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    // Member Functions

        //- Turbulent dispersion force (includes drag already)
        //virtual tmp<volVectorField> F() const = 0;

	// Diffusivity
        virtual tmp<volScalarField> D() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
