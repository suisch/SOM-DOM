
{

	Info << "Reading solution species dictionary" << endl;
	const dictionary& solutionSpeciesDict = twoPhaseProperties.subDict("solutionSpecies");
	Info << "Reading solution species name\n" << endl;
	speciesTable solutionSpecies
	(
		solutionSpeciesDict.toc()
	);

	scalarList maxY(solutionSpecies.size());
	forAll(solutionSpecies,i)
	{
		maxY[i] = readScalar(mesh.solutionDict().subDict("PISO").subDict("maxConcentration").lookup(solutionSpecies[i]));
	}

	Info << "Reading compressive flux coefficient\n" << endl;	
	scalar cYi(readScalar(mesh.solutionDict().subDict("PISO").lookup("cYi")));

	//Transport step
	Info << "transport step" << endl;
	word alpharScheme("div(phirb,alpha)");
	word YiScheme("div(phi,Yi)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(cYi*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    surfaceScalarField fluxDir = fvc::snGrad(alpha1)*mesh.magSf();
	surfaceScalarField alphaUp = upwind<scalar>(mesh,fluxDir).interpolate(alpha1);
	surfaceScalarField alphaDown = downwind<scalar>(mesh,fluxDir).interpolate(alpha1);

	{//C_glu transport
		volScalarField& Yi = C_glu;
		dimensionedScalar D1Yi = D_glu_Y1;
		dimensionedScalar D2Yi = D_glu_Y2;
		dimensionedScalar HYi = H_glu;

		Yi.correctBoundaryConditions();

		//compute diffusion coefficient
		surfaceScalarField D = fvc::interpolate(D1Yi*alpha1+HYi*D2Yi*(1-alpha1))/fvc::interpolate(alpha1+HYi*(1-alpha1));

		surfaceScalarField phiH = -D*(1-HYi)/fvc::interpolate((alpha1+HYi*(1-alpha1)))*fvc::snGrad(alpha1)*mesh.magSf();

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

		surfaceScalarField YiUp   = upwind<scalar>(mesh,fluxDir).interpolate(Yi);
		surfaceScalarField YiDown = downwind<scalar>(mesh,fluxDir).interpolate(Yi);
        surfaceScalarField deltaYi1 = max(-max(Yi),min(max(Yi),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
		surfaceScalarField deltaYi2 = YiDown/(alphaDown+(1-alphaDown)*HYi)-HYi*YiUp/(alphaUp+(1-alphaUp)*HYi);

		//interface compression flux
		surfaceScalarField phiY =
        fvc::flux
        (
            phi,
            Yi,
            YiScheme
        )
		+phiComp*sign(deltaYi2)*max(mag(deltaYi1),mag(deltaYi2));

    	MULES::explicitSolve(Yi, phi, phiY, maxY[0], 0);

		volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

		while (simple.correctNonOrthogonal())
		{
			//compute Y equation
			fvScalarMatrix YiEqn
			(
				fvm::ddt(Yi)
		        - fvm::laplacian(D,Yi)
				- fvm::div(phiH,Yi,"div(phiH,Yi)")
				==
				SpY
			);

			//solve equations
		    if (simple.finalNonOrthogonalIter())
		    {
		        YiEqn.solve(mesh.solutionDict().solver("YiFinal"));
		    }
		    else
		    {
				YiEqn.solve(mesh.solutionDict().solver("Yi"));
			}
		}

		Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
		<< "  Min(Yi) = " << gMin(Yi.internalField())
		<< "  Max(Yi) = " << gMax(Yi.internalField())
		<< endl;
	}
	{//C_DOM transport
		volScalarField& Yi = C_DOM;
		dimensionedScalar D1Yi = D_DOM_Y1;
		dimensionedScalar D2Yi = D_DOM_Y2;
		dimensionedScalar HYi = H_DOM;

		Yi.correctBoundaryConditions();

		//compute diffusion coefficient
		surfaceScalarField D = fvc::interpolate(D1Yi*alpha1+HYi*D2Yi*(1-alpha1))/fvc::interpolate(alpha1+HYi*(1-alpha1));

		surfaceScalarField phiH = -D*(1-HYi)/fvc::interpolate((alpha1+HYi*(1-alpha1)))*fvc::snGrad(alpha1)*mesh.magSf();

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

		surfaceScalarField YiUp   = upwind<scalar>(mesh,fluxDir).interpolate(Yi);
		surfaceScalarField YiDown = downwind<scalar>(mesh,fluxDir).interpolate(Yi);
        surfaceScalarField deltaYi1 = max(-max(Yi),min(max(Yi),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
		surfaceScalarField deltaYi2 = YiDown/(alphaDown+(1-alphaDown)*HYi)-HYi*YiUp/(alphaUp+(1-alphaUp)*HYi);

		//interface compression flux
		surfaceScalarField phiY =
        fvc::flux
        (
            phi,
            Yi,
            YiScheme
        )
		+phiComp*sign(deltaYi2)*max(mag(deltaYi1),mag(deltaYi2));

    	MULES::explicitSolve(Yi, phi, phiY, maxY[1], 0);

		volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

		while (simple.correctNonOrthogonal())
		{
			//compute Y equation
			fvScalarMatrix YiEqn
			(
				fvm::ddt(Yi)
		        - fvm::laplacian(D,Yi)
				- fvm::div(phiH,Yi,"div(phiH,Yi)")
				==
				SpY
			);

			//solve equations
		    if (simple.finalNonOrthogonalIter())
		    {
		        YiEqn.solve(mesh.solutionDict().solver("YiFinal"));
		    }
		    else
		    {
				YiEqn.solve(mesh.solutionDict().solver("Yi"));
			}
		}

		Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
		<< "  Min(Yi) = " << gMin(Yi.internalField())
		<< "  Max(Yi) = " << gMax(Yi.internalField())
		<< endl;
	}
	{//C_O2 transport
		volScalarField& Yi = C_O2;
		dimensionedScalar D1Yi = D_O2_Y1;
		dimensionedScalar D2Yi = D_O2_Y2;
		dimensionedScalar HYi = H_O2;

		Yi.correctBoundaryConditions();

		//compute diffusion coefficient
		surfaceScalarField D = fvc::interpolate(D1Yi*alpha1+HYi*D2Yi*(1-alpha1))/fvc::interpolate(alpha1+HYi*(1-alpha1));

		surfaceScalarField phiH = -D*(1-HYi)/fvc::interpolate((alpha1+HYi*(1-alpha1)))*fvc::snGrad(alpha1)*mesh.magSf();

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

		surfaceScalarField YiUp   = upwind<scalar>(mesh,fluxDir).interpolate(Yi);
		surfaceScalarField YiDown = downwind<scalar>(mesh,fluxDir).interpolate(Yi);
        surfaceScalarField deltaYi1 = max(-max(Yi),min(max(Yi),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
		surfaceScalarField deltaYi2 = YiDown/(alphaDown+(1-alphaDown)*HYi)-HYi*YiUp/(alphaUp+(1-alphaUp)*HYi);

		//interface compression flux
		surfaceScalarField phiY =
        fvc::flux
        (
            phi,
            Yi,
            YiScheme
        )
		+phiComp*sign(deltaYi2)*max(mag(deltaYi1),mag(deltaYi2));

    	MULES::explicitSolve(Yi, phi, phiY, maxY[2], 0);

		volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

		while (simple.correctNonOrthogonal())
		{
			//compute Y equation
			fvScalarMatrix YiEqn
			(
				fvm::ddt(Yi)
		        - fvm::laplacian(D,Yi)
				- fvm::div(phiH,Yi,"div(phiH,Yi)")
				==
				SpY
			);

			//solve equations
		    if (simple.finalNonOrthogonalIter())
		    {
		        YiEqn.solve(mesh.solutionDict().solver("YiFinal"));
		    }
		    else
		    {
				YiEqn.solve(mesh.solutionDict().solver("Yi"));
			}
		}

		Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
		<< "  Min(Yi) = " << gMin(Yi.internalField())
		<< "  Max(Yi) = " << gMax(Yi.internalField())
		<< endl;
	}
	{//C_CO2 transport
		volScalarField& Yi = C_CO2;
		dimensionedScalar D1Yi = D_CO2_Y1;
		dimensionedScalar D2Yi = D_CO2_Y2;
		dimensionedScalar HYi = H_CO2;

		Yi.correctBoundaryConditions();

		//compute diffusion coefficient
		surfaceScalarField D = fvc::interpolate(D1Yi*alpha1+HYi*D2Yi*(1-alpha1))/fvc::interpolate(alpha1+HYi*(1-alpha1));

		surfaceScalarField phiH = -D*(1-HYi)/fvc::interpolate((alpha1+HYi*(1-alpha1)))*fvc::snGrad(alpha1)*mesh.magSf();

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

		surfaceScalarField YiUp   = upwind<scalar>(mesh,fluxDir).interpolate(Yi);
		surfaceScalarField YiDown = downwind<scalar>(mesh,fluxDir).interpolate(Yi);
        surfaceScalarField deltaYi1 = max(-max(Yi),min(max(Yi),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
		surfaceScalarField deltaYi2 = YiDown/(alphaDown+(1-alphaDown)*HYi)-HYi*YiUp/(alphaUp+(1-alphaUp)*HYi);

		//interface compression flux
		surfaceScalarField phiY =
        fvc::flux
        (
            phi,
            Yi,
            YiScheme
        )
		+phiComp*sign(deltaYi2)*max(mag(deltaYi1),mag(deltaYi2));

    	MULES::explicitSolve(Yi, phi, phiY, maxY[3], 0);

		volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

		while (simple.correctNonOrthogonal())
		{
			//compute Y equation
			fvScalarMatrix YiEqn
			(
				fvm::ddt(Yi)
		        - fvm::laplacian(D,Yi)
				- fvm::div(phiH,Yi,"div(phiH,Yi)")
				==
				SpY
			);

			//solve equations
		    if (simple.finalNonOrthogonalIter())
		    {
		        YiEqn.solve(mesh.solutionDict().solver("YiFinal"));
		    }
		    else
		    {
				YiEqn.solve(mesh.solutionDict().solver("Yi"));
			}
		}

		Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
		<< "  Min(Yi) = " << gMin(Yi.internalField())
		<< "  Max(Yi) = " << gMax(Yi.internalField())
		<< endl;
	}
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // fvScalarMatrix C_gluEqn
    // (
	// fvm::ddt(C_glu)
	// + fvm::div(phi, C_glu, "div(phi,C)")
	// - fvm::laplacian(D_glu, C_glu, "laplacian(C)")
    // );

    // C_gluEqn.solve();

    // fvScalarMatrix C_DOMEqn
    // (
	// fvm::ddt(C_DOM)
	// + fvm::div(phi, C_DOM, "div(phi,C)")
	// - fvm::laplacian(D_DOM, C_DOM, "laplacian(C)")
    // );

    // C_DOMEqn.solve();

	// fvScalarMatrix C_O2Eqn
    // (
	// fvm::ddt(C_O2)
	// + fvm::div(phi, C_O2, "div(phi,C)")
	// - fvm::laplacian(D_O2, C_O2, "laplacian(C)")
    // );

    // C_O2Eqn.solve();

    // fvScalarMatrix C_CO2Eqn
    // (
	// fvm::ddt(C_CO2)
	// + fvm::div(phi, C_CO2, "div(phi,C)")
	// - fvm::laplacian(D_CO2, C_CO2, "laplacian(C)")
    // );

    // C_CO2Eqn.solve();

    //  fvScalarMatrix C_mAEEqn
    // (
	// fvm::ddt(C_mAE)
	// + fvm::div(phi, C_mAE, "div(phi,C)")
	// - fvm::laplacian(DT, C_mAE, "laplacian(C)")
    // );

    // C_mAEEqn.solve();
	
	// fvScalarMatrix C_mF1Eqn
    // (
	// fvm::ddt(C_mF1)
	// + fvm::div(phi, C_mF1, "div(phi,C)")
	// - fvm::laplacian(DT, C_mF1, "laplacian(C)")
    // );

    // C_mF1Eqn.solve();

	// fvScalarMatrix C_mF2Eqn
    // (
	// fvm::ddt(C_mF2)
	// + fvm::div(phi, C_mF2, "div(phi,C)")
	// - fvm::laplacian(DT, C_mF2, "laplacian(C)")
    // );

    // C_mF2Eqn.solve();
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//reaction step
	Info << "Reaction Step" << endl;

	double celli = 0;
	
	forAll(C_mAE.boundaryField(), patchi)
	{

		if (Surf.boundaryField()[patchi].type()=="reactingWall")
		{

			const labelList& cellOwner = C_mAE.boundaryField()[patchi].patch().faceCells();
			
			forAll(C_mAE.boundaryField()[patchi], facei)
			{
				
				celli = cellOwner[facei];
				//Equilibrium reaction for SOM-DOM
				double C_DOM_temp = C_DOM[celli];
				double C_SOM_temp = C_SOM[celli];

				C_DOM[celli] = kD_SOMDOM.value()/(1+kD_SOMDOM.value())*(C_DOM_temp + C_SOM_temp);
				C_SOM[celli] = 1/(1+kD_SOMDOM.value())*(C_DOM_temp + C_SOM_temp);

				C_DOM_temp = NULL;
				C_SOM_temp = NULL;

				double time_step = runTime.deltaT().value();
				
				//set concentration for BRNS module
				//for (int i = 0 ; i < numComp; i++){
				concentrationTemp_[0] = C_glu[celli];
				concentrationTemp_[1] = C_DOM[celli];
				concentrationTemp_[2] = 0;				//C_O2
				concentrationTemp_[3] = 0;				//C_CO2
				concentrationTemp_[4] = C_mAE[celli];
				concentrationTemp_[5] = C_mF1[celli];
				concentrationTemp_[6] = C_mF2[celli];

				for (int i = 0 ; i < numComp; i++)
				{
					boundary_flag[i] = 0;
				}

				/*cout << "concentration A = " << concentrationTemp_[0] << endl;
				cout << "concentration B = " << concentrationTemp_[1] << endl;
				cout << "concentration C = " << concentrationTemp_[2] << endl;
				cout << "concentration bio = " << concentrationTemp_[3] << endl;*/

				// Run BRNS;
				invokebrns_
				(
					concentrationTemp_,
					concentrationTemp_,
					concentrationTemp_,
					&numComp,
					&time_step,
					boundary_flag,
					&return_value,
					&pos_x,
					&pos_y,
					&pos_z,
					&porosity,
					&waterSaturation,
					NULL
				);

				if (concentrationTemp_[2] > 1e-20)
				{ 
					C_O2[celli] += concentrationTemp_[2];
				}
				else if ( C_O2[celli] > 1e-20)
				{
					C_O2[celli] += 0;
				}
				else
				{
					C_O2[celli] = 0;
				}

				if (concentrationTemp_[3] > 1e-20)
				{ 
					C_CO2[celli] += concentrationTemp_[3];
				}
				else if ( C_CO2[celli] > 1e-20)
				{
					C_CO2[celli] += 0;
				}
				else
				{
					C_CO2[celli] = 0;
				}

				C_glu[celli]    = concentrationTemp_[0];//mol/m3
				C_DOM[celli]    = concentrationTemp_[1];//mol/m3
				// C_O2[celli]     += concentrationTemp_[2];//mol/m3
				// C_CO2[celli]    += concentrationTemp_[3];//mol/m3
				C_mAE[celli]    = concentrationTemp_[4];//mol/m3
				C_mF1[celli]    = concentrationTemp_[5];//mol/m3
				C_mF2[celli]    = concentrationTemp_[6];//mol/m3

			}
		}
	}		
}

Info << "End of reaction update." << endl;
Info << endl;