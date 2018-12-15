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
    Foam::multiphaseTurbulenceModel

Description

SourceFiles
    multiphaseTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseTurbulenceModel_H
#define multiphaseTurbulenceModel_H

#include "volFields.H"
#include "dragModel.H"
#include "phaseModel.H"
#include "nearWallDist.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiphaseTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class multiphaseTurbulenceModel
{

protected:

    // Protected but accessible 

        const phaseModel& phase1_;
        const phaseModel& phase2_;
        const volVectorField& U1_;
        const volVectorField& U2_;
        const volScalarField& alpha1_;
        const volScalarField& alpha2_;
        const surfaceScalarField& phi1_;
        const surfaceScalarField& phi2_;
        const dragModel& drag1_;
        const dimensionedScalar& rho1_;
        const dimensionedScalar& rho2_;
        const dimensionedScalar& da_;
        const dimensionedScalar& nu1_;
        const dimensionedScalar& nu2_;

        nearWallDist y_;

private:

        //- Disallow default bitwise copy construct
        multiphaseTurbulenceModel(const multiphaseTurbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const multiphaseTurbulenceModel&);


public:

    //- Runtime type information
    TypeName("multiphaseTurbulenceModel");

    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            multiphaseTurbulenceModel,
            multiphaseTurbulenceModel,
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
        multiphaseTurbulenceModel
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dragModel& drag1
        );

    // Selectors
    
        //- Return ref to turbulence model
        static autoPtr<multiphaseTurbulenceModel> New
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dragModel& drag1
        );

    //- Destructor
    virtual ~multiphaseTurbulenceModel();

    // Member Functions

        virtual void solve() = 0;

        // Access to turbulence
        virtual tmp<volScalarField> k() const = 0;

        // Access to dissipation 
        virtual tmp<volScalarField> epsilon() const = 0;

        // Access to turbulence viscosity
        virtual tmp<volScalarField> nut() const = 0;

        // Access to effective viscosity // may implement this 
        //virtual tmp<volScalarField> nuEff() const = 0;

        // read dict for consts
        //virtual bool read() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
