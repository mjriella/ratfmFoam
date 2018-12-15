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
    Foam::dispersedv2f

Description
	Particle phase v2f model. Wall boundary conditions are:

		k       = zeroGradient/fixedValue
        epsilon = hard coded
        v2      = zeroGradient/fixedValue
        f       = fixedValue 

Note
	The model is sensitive to inital conditions and the BC depends on flow
	regime.

References
	Riella, M., Tabor, G., and Kahraman, R. (2018b). Inhomogeneity and 
	anisotropy in Eulerian-	Eulerian near-wall modelling. International 
	Journal of Multiphase Flow.

SourceFiles
    dispersedv2f.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedv2f_H
#define dispersedv2f_H

#include "particleTurbulence.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseTurbulenceModels
{

/*---------------------------------------------------------------------------*\
                           Class dispersedv2f Declaration
\*---------------------------------------------------------------------------*/

class dispersedv2f
:
    public particleTurbulence
{

protected:

        dimensionedScalar Cmu_;
        dimensionedScalar CmuKEps_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar CL_;
        dimensionedScalar Ceta_;
        dimensionedScalar Ceps2_;
        dimensionedScalar sigmaK_;
        dimensionedScalar sigmaEps_;

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField nut_;

		// wall normal component 
		volScalarField v1_;

		// damping function
		volScalarField f_;

        dimensionedScalar v1Min_;
        dimensionedScalar fMin_;

    // Protected Member Functions

        //- Apply Davidson correction to nut
        tmp<volScalarField> davidsonCorrectNut
        (
            const tmp<volScalarField>& value
        ) const;

        //- Return time scale, Ts
        tmp<volScalarField> Ts() const;

        //- Return length scale, Ls
        tmp<volScalarField> Ls() const;

public:

    //- Runtime type information
    TypeName("dispersedv2f");


    // Constructors

        //- Construct from components
        dispersedv2f
		(
			//const dictionary& dict,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dragModel& drag1
		);


    //- Destructor
    virtual ~dispersedv2f();


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

        // call wall normal in phase 2
        virtual tmp<volScalarField> v2() const
        {
            return v1_;
        }

        // call wall damping in phase 2
        virtual tmp<volScalarField> f() const
        {
            return f_;
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
