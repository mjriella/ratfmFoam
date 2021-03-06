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
    Foam::fluidTurbulence

Description
	Base class for fluid phase turbulence models.

SourceFiles
    fluidTurbulence.C

\*---------------------------------------------------------------------------*/

#ifndef fluidTurbulence_H
#define fluidTurbulence_H

#include "../multiphaseTurbulenceModel/multiphaseTurbulenceModel.H"
#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseTurbulenceModels
{

/*---------------------------------------------------------------------------*\
                           Class fluidTurbulence Declaration
\*---------------------------------------------------------------------------*/

class fluidTurbulence
:   
    public multiphaseTurbulenceModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        fluidTurbulence(const fluidTurbulence&);

        //- Disallow default bitwise assignment
        void operator=(const fluidTurbulence&);


protected:

    // Protected data

        IOdictionary multiphaseTurbulenceProperties_;

        dictionary coeffsDict_;

        // limiters
        dimensionedScalar kMin_;
        dimensionedScalar epsilonMin_;

public:

    //- Runtime type information
    TypeName("fluidTurbulence");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidTurbulence,
        dictionary,
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dragModel& drag1
        ),
        (phase1, phase2, alpha1, alpha2, drag1)
    );


    // Constructors

        //- Construct from components
        fluidTurbulence
        (
            //const dictionary& dict,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dragModel& drag1
        );


    // Selectors

        static autoPtr<fluidTurbulence> New
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dragModel& drag1
        );


    //- Destructor
    virtual ~fluidTurbulence();


    // Member Functions

        // Access
        const dimensionedScalar kMin() const
        {
            return kMin_;
        } 

        const dimensionedScalar epsilonMin() const
        {
            return epsilonMin_;
        }

        dimensionedScalar kMin() 
        {
            return kMin_;
        } 

        dimensionedScalar epsilonMin() 
        {
            return epsilonMin_;
        }

        // access to dict
        virtual const dictionary& coeffsDict() const
        {
            return coeffsDict_;
        }

        virtual void solve();
        
        // Read the multiphase dict
        //virtual void read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multiphaseTurbulenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
