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
    Foam::kEpsilon

Description
	Fluid phase kEpsilon model with co-variance coupling. 

References
	Riella, M., Kahraman, R., and Tabor, G. R. (2018). Reynolds-averaged two-
	fluid model prediction of moderately dilute gas-solid flow over a backward
	-facing step. International Journal of Multiphase Flow, 106:95 – 108.

    Fox, R. O (2014). On multiphase turbulence models for collisional fluid-
    particle flows. Journal of Fluid Mechanics, 742:368-424.

SourceFiles
    kEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilon_H
#define kEpsilon_H

#include "fluidTurbulence.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseTurbulenceModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilon Declaration
\*---------------------------------------------------------------------------*/

class kEpsilon
:
    public fluidTurbulence
{

protected:

        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar sigmaEps_;
        dimensionedScalar kappa_;
        dimensionedScalar E_;

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField nut_;

    // Protected Member Functions

		// Enforce wall func
		void correctEpsilon();
		
		// Enforce nut wall func
		void correctNut();

public:

    //- Runtime type information
    TypeName("kEpsilon");


    // Constructors

        //- Construct from components
        kEpsilon
		(
			//const dictionary& dict,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dragModel& drag1
		);


    //- Destructor
    virtual ~kEpsilon();


    // Member functions
        
        // solve syst of eqs
        virtual void solve();

        virtual bool read();

        // call turbulent kinetic energy in phase 2
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        // call dissipation in phase 2
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }
        // call turbulence viscosity in phase 2
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multiphaseTurbulenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
