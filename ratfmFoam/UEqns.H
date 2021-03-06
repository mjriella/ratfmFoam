fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

// Declare turbulence viscosity
volScalarField nut1 = multiphaseTurbulence1->nut();
volScalarField nut2 = multiphaseTurbulence2->nut();
volScalarField nuEff1 = nut1;
volScalarField nuEff2 = nut2 + nu2;

{
    {
        volTensorField gradU1T(T(fvc::grad(U1)));

        if (kineticTheory.on())
        {
            kineticTheory.solve(gradU1T);
            nuEff1 += kineticTheory.mu1()/rho1;
        }
        else
        {
            nuEff1 = nut1;
        }

        // - Particle stress tensor
        volTensorField Rc1
        (
            "Rc1",
            (((2.0/3.0)*I)*nuEff1)*tr(gradU1T) - nuEff1*gradU1T
        );

        surfaceScalarField phiR1
        (
            -fvc::interpolate(nuEff1)*mesh.magSf()*fvc::snGrad(alpha1)
            /fvc::interpolate(alpha1 + alphaSmall)
        );

        U1Eqn =
        (
            fvm::ddt(U1)
          + fvm::div(phi1, U1, "div(phi1,U1)")
          - fvm::Sp(fvc::div(phi1), U1)
		  - fvm::laplacian(nuEff1, U1) 
          + fvc::div(Rc1)
          + fvm::div(phiR1, U1, "div(phi1,U1)") 
          - fvm::Sp(fvc::div(phiR1), U1)
          + (fvc::grad(alpha1)/(fvc::average(alpha1) + alphaSmall) & Rc1)
         ==
          - fvm::Sp(alpha2/rho1*K, U1)
        );
        U1Eqn.relax();
    }

    {
        volTensorField gradU2T(T(fvc::grad(U2)));

        // - Fluid phase stress tensor
        volTensorField Rc2
        (
            "Rc2",
            (((2.0/3.0)*I)*nuEff2)*tr(gradU2T) - nuEff2*gradU2T
        );

        surfaceScalarField phiR2
        (
            -fvc::interpolate(nuEff2)*mesh.magSf()*fvc::snGrad(alpha2)
            /fvc::interpolate(alpha2 + alphaSmall)
        );

        U2Eqn =
        (
            fvm::ddt(U2)
          + fvm::div(phi2, U2, "div(phi2,U2)")
          - fvm::Sp(fvc::div(phi2), U2)
          - fvm::laplacian(nuEff2, U2) 
          + fvc::div(Rc2)
          + fvm::div(phiR2, U2, "div(phi2,U2)") 
          - fvm::Sp(fvc::div(phiR2), U2)
          + (fvc::grad(alpha2)/(fvc::average(alpha2) + alphaSmall) & Rc2)
         ==
          - fvm::Sp(alpha1/rho2*K, U2)
        );
        U2Eqn.relax();
    }
}
