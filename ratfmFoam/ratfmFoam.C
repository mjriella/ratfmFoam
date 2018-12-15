/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    ratfmFoam 

Description
    A Reynolds-Averaged Two-Fluid model for turbulent fluid-particle flows.
	The solver features a multiscale approach to account for the partitioning
	of particle energy. Multiphase turbulence modelling is accounted for and 
	enables coupling across phases through the co-variance term. 

Note 
	Particle phase = 1 and fluid phase = 2 with turbulence always on. 
	Additionally, the turbulence models should match for consistency and
	their boundary conditions are hard coded.

Author
    Matthew Riella, University of Exeter
    
References
	Riella, M., Kahraman, R., and Tabor, G. R. (2018). Reynolds-averaged two-
	fluid model prediction of moderately dilute gas-solid flow over a backward
	-facing step. International Journal of Multiphase Flow, 106:95 – 108.

    Fox, R. O (2014). On multiphase turbulence models for collisional fluid-
    particle flows. Journal of Fluid Mechanics, 742:368-424.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "fixedValueFvsPatchFields.H"
#include "turbulentDispersionModel.H"
#include "dragModel.H"
#include "phaseModel.H"
#include "kineticTheoryModel.H"
#include "multiphaseTurbulenceModel.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "readPPProperties.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readratfmFoamControls.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaEqn.H"
            #include "dragCoeffs.H"
            #include "UEqns.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
				multiphaseTurbulence1->solve();
				multiphaseTurbulence2->solve();
            }
        }

        #include "write.H"
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
