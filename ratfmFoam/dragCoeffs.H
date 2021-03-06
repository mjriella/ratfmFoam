    volVectorField Ur(U1 - U2);
    volScalarField magUr(mag(Ur));

    volScalarField Ka(drag1->K(magUr));
    volScalarField K(Ka);

    // bit of a cheat for Zuber dispersion model
	volScalarField beta
	(
    	IOobject
    	(
        	"beta",
        	runTime.timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::NO_WRITE
    	),
    	alpha1*(1.0-alpha1)*K
	);

    // Remove lift and drag at fixed-flux boundaries
    forAll(phi1.boundaryField(), patchi)
    {
        if (isA<fixedValueFvsPatchScalarField>(phi1.boundaryField()[patchi]))
        {
            K.boundaryField()[patchi] = 0.0;
        }
    }

