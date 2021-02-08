{

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

	//Re-initialize transfer flux
	phiD *=0;

    forAll(solutionSpecies, i)
	{
		volScalarField& Yi = Y[i];
		dimensionedScalar D1Yi = D1Y[i];
		dimensionedScalar D2Yi = D2Y[i];
		dimensionedScalar HYi = HY[i];

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

    	MULES::explicitSolve(Yi, phi, phiY, maxY[i], 0);

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

		phiD+=Mw[i]*
			 (
			 	D*fvc::snGrad(Yi)*mesh.magSf()
			   +fvc::flux(phiH,Yi,"div(phiH,Yi)")
			 );
	}

	Mflux = fvc::div(phiD*alphaDown)-alpha1*fvc::div(phiD);

	//compute Yi1 and Yi2
	forAll(solutionSpecies, i)
	{
		volScalarField& Yi = Y[i];
		volScalarField& Y1i = Y1[i];
		volScalarField& Y2i = Y2[i];
		dimensionedScalar HYi = HY[i];
		Y1i = Yi/(alpha1+HYi*(1-alpha1));
		Y2i = HYi*Yi/(alpha1+HYi*(1-alpha1));
	}

	forAll(solutionSpecies, i)
	{
		volScalarField& Yi = Y[i];
		volScalarField& Y1i = Y1[i];
		volScalarField& Y2i = Y2[i];
		Yi = Y1i*alpha1 + Y2i*(1-alpha1);
	}

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	//reaction step
	Info << "Reaction Step" << endl;

	double celli = 0;
	
	forAll(Y[0].boundaryField(), patchi)
	{
		
		if (Surf.boundaryField()[patchi].type()=="reactingWall")
		{

			const labelList& cellOwner = Y[0].boundaryField()[patchi].patch().faceCells();
			
			forAll(Y[0].boundaryField()[patchi], facei)
			{
				
				celli = cellOwner[facei];
				//Equilibrium reaction for SOM-DOM
				{
					volScalarField& C_DOM_temp = Y[1]; //C_DOM
					double C_DOM_temp_celli = C_DOM_temp[celli];
					double C_SOM_temp = C_SOM[celli];

					C_DOM_temp[celli] = kD_SOMDOM.value()/(1+kD_SOMDOM.value())*(C_DOM_temp_celli + C_SOM_temp);
					C_SOM[celli] = 1/(1+kD_SOMDOM.value())*(C_DOM_temp_celli + C_SOM_temp);
				}

				double time_step = runTime.deltaT().value();
				
				//set concentration for BRNS module
				//for (int i = 0 ; i < numComp; i++){
				// int iCounter = 0;
				// forAll(solutionSpecies, i)
				// {
				// 	concentrationTemp_[iCounter] = Y[i][celli];
				// 	iCounter += 1;
				// }
				// iCounter = NULL;

				//concentrationTemp_[0] = C_glu[celli];
				//concentrationTemp_[1] = C_DOM[celli];
				//concentrationTemp_[2] = 0;				//C_O2
				//concentrationTemp_[3] = 0;				//C_CO2
				concentrationTemp_[0] = Y[0][celli];
				concentrationTemp_[1] = Y[1][celli];
				concentrationTemp_[2] = 0;				//C_O2
				concentrationTemp_[3] = 0;				//C_CO2
				concentrationTemp_[4] = C_mAE[celli];
				concentrationTemp_[5] = C_mF1[celli];
				concentrationTemp_[6] = C_mF2[celli];

				for (int i = 0 ; i < numComp; i++)
				{
					boundary_flag[i] = 0;
				}

				/*boundary_flag[0] = 0;
				boundary_flag[1] = 0;
				boundary_flag[2] = 0;
				boundary_flag[3] = 0;
				boundary_flag[4] = 0;
				boundary_flag[5] = 0;
				boundary_flag[6] = 0;*/

				/*cout << "concentration A = " << concentrationTemp_[0] << endl;
				cout << "concentration B = " << concentrationTemp_[1] << endl;
				cout << "concentration C = " << concentrationTemp_[2] << endl;
				cout << "concentration bio = " << concentrationTemp_[3] << endl;*/

		// Info << "Invoking BRNS" << endl;
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

				// cout << "concentration C_glu = " << concentrationTemp_[0] << endl;
				// cout << "concentration C_DOM = " << concentrationTemp_[1] << endl;
				// cout << "concentration C_O2 = " << concentrationTemp_[2] << endl;
				// cout << "concentration C_CO2 = " << concentrationTemp_[3] << endl;
				// cout << "concentration C_mAE = " << concentrationTemp_[4] << endl;
				// cout << "concentration C_mF1 = " << concentrationTemp_[5] << endl;
				// cout << "concentration C_mF2 = " << concentrationTemp_[6] << endl;

				if (concentrationTemp_[2] > 1e-20)
				{ 
					Y[2][celli] += concentrationTemp_[2];
				}
				else if ( Y[2][celli] > 1e-20)
				{
					Y[2][celli] += 0;
				}
				else
				{
					Y[2][celli] = 0;
				}

				if (concentrationTemp_[3] > 1e-20)
				{ 
					Y[3][celli] += concentrationTemp_[3];
				}
				else if ( Y[3][celli] > 1e-20)
				{
					Y[3][celli] += 0;
				}
				else
				{
					Y[3][celli] = 0;
				}

				if (concentrationTemp_[4] > 1e-20)
				{ 
					C_mAE[celli] += concentrationTemp_[4];
				}
				else if ( C_mAE[celli] > 1e-20)
				{
					C_mAE[celli] += 0;
				}
				else
				{
					C_mAE[celli] = 0;
				}

				if (concentrationTemp_[5] > 1e-20)
				{ 
					C_mF1[celli] += concentrationTemp_[5];
				}
				else if ( C_mF1[celli] > 1e-20)
				{
					C_mF1[celli] += 0;
				}
				else
				{
					C_mF1[celli] = 0;
				}

				if (concentrationTemp_[6] > 1e-20)
				{ 
					C_mF2[celli] += concentrationTemp_[6];
				}
				else if ( C_mF2[celli] > 1e-20)
				{
					C_mF2[celli] += 0;
				}
				else
				{
					C_mF2[celli] = 0;
				}

				Y[0][celli]    = concentrationTemp_[0];//mol/m3
				Y[1][celli]    = concentrationTemp_[1];//mol/m3
				// Y[2][celli]    = concentrationTemp_[2];//mol/m3
				// Y[3][celli]    = concentrationTemp_[3];//mol/m3
				//C_glu[celli]    = concentrationTemp_[0];//mol/m3
				//C_DOM[celli]    = concentrationTemp_[1];//mol/m3
				//C_O2[celli]     += concentrationTemp_[2];//mol/m3
				//C_CO2[celli]    += concentrationTemp_[3];//mol/m3
				// C_mAE[celli]    = concentrationTemp_[4];//mol/m3
				// C_mF1[celli]    = concentrationTemp_[5];//mol/m3
				// C_mF2[celli]    = concentrationTemp_[6];//mol/m3
				// cout << "concentration C_glu = " << Y[0][celli] << endl;
				// cout << "concentration C_DOM = " << Y[1][celli] << endl;
				// cout << "concentration C_O2 = " << Y[2][celli] << endl;
				// cout << "concentration C_CO2 = " << Y[3][celli] << endl;
				// cout << "concentration C_mAE = " << C_mAE[celli] << endl;
				// cout << "concentration C_mF1 = " << C_mF1[celli] << endl;
				// cout << "concentration C_mF2 = " << C_mF2[celli] << endl;
			}
		}
	}		
}

Info << "End of reaction step" << endl;
Info << " " << endl;
