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
    Foam::kineticTheoryModel

Description

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
{
    // Private data

        const phaseModel& phase1_;
        const volVectorField& U1_;
        const volVectorField& U2_;
        const volScalarField& alpha1_;
        const surfaceScalarField& phi1_;

        const dragModel& drag1_;

        const dimensionedScalar& rho1_;
        const dimensionedScalar& da_;
        const dimensionedScalar& nu1_;

        //- dictionary holding the modeling info
        IOdictionary kineticTheoryProperties_;

        //- use kinetic theory or not.
        Switch kineticTheory_;

        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

        autoPtr<conductivityModel> conductivityModel_;

        autoPtr<kineticTheoryModels::radialModel> radialModel_;

        autoPtr<granularPressureModel> granularPressureModel_;

        //- coefficient of restitution
        const dimensionedScalar e_;

        //- maximum packing
        const dimensionedScalar alphaMax_;

        //- Uncorrelated energy
        volScalarField Theta_;

        //- The granular viscosity
        volScalarField mu1_;

        //- The granular viscosity
        volScalarField nuKT_;

        //- The granular pressure
        volScalarField pa_;

        //- The granular temperature conductivity
        volScalarField kappa_;

        //- The radial distribution function
        volScalarField gs0_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);


public:

    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const phaseModel& phase1,
            const volVectorField& U2,
            const volScalarField& alpha1,
            const dragModel& drag1
        );


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions

        void solve(const volTensorField& gradU1t);

        bool on() const
        {
            return kineticTheory_;
        }

        const volScalarField& mu1() const
        {
            return mu1_;
        }

        const volScalarField& nuKT() const
        {
            return nuKT_;
        }

        const volScalarField& pa() const
        {
            return pa_;
        }

        const volScalarField& kappa() const
        {
            return kappa_;
        }
        
        const volScalarField& Theta() const
        {
            return Theta_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
