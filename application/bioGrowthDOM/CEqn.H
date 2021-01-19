

{
    fvScalarMatrix C_gluEqn
    (
	fvm::ddt(C_glu)
	+ fvm::div(phi, C_glu, "div(phi,C)")
	- fvm::laplacian(D_glu, C_glu, "laplacian(C)")
    );

    C_gluEqn.solve();

    fvScalarMatrix C_DOMEqn
    (
	fvm::ddt(C_DOM)
	+ fvm::div(phi, C_DOM, "div(phi,C)")
	- fvm::laplacian(D_DOM, C_DOM, "laplacian(C)")
    );

    C_DOMEqn.solve();

	fvScalarMatrix C_O2Eqn
    (
	fvm::ddt(C_O2)
	+ fvm::div(phi, C_O2, "div(phi,C)")
	- fvm::laplacian(D_O2, C_O2, "laplacian(C)")
    );

    C_O2Eqn.solve();

    fvScalarMatrix C_CO2Eqn
    (
	fvm::ddt(C_CO2)
	+ fvm::div(phi, C_CO2, "div(phi,C)")
	- fvm::laplacian(D_CO2, C_CO2, "laplacian(C)")
    );

    C_CO2Eqn.solve();

/*    fvScalarMatrix C_mAEEqn
    (
	fvm::ddt(C_mAE)
	+ fvm::div(phi, C_mAE, "div(phi,C)")
	- fvm::laplacian(DT, C_mAE, "laplacian(C)")
    );

    C_mAEEqn.solve();
	
	fvScalarMatrix C_mF1Eqn
    (
	fvm::ddt(C_mF1)
	+ fvm::div(phi, C_mF1, "div(phi,C)")
	- fvm::laplacian(DT, C_mF1, "laplacian(C)")
    );

    C_mF1Eqn.solve();

	fvScalarMatrix C_mF2Eqn
    (
	fvm::ddt(C_mF2)
	+ fvm::div(phi, C_mF2, "div(phi,C)")
	- fvm::laplacian(DT, C_mF2, "laplacian(C)")
    );

    C_mF2Eqn.solve();*/
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	double celli = 0;
	
	forAll(C_glu.boundaryField(), patchi)
	{
		if (Surf.boundaryField()[patchi].type()=="reactingWall")
		{	
			const labelList& cellOwner = C_glu.boundaryField()[patchi].patch().faceCells();
			
			forAll(C_glu.boundaryField()[patchi], facei)
			{
				//Equilibrium reaction for SOM-DOM
				double C_DOM_temp = C_DOM[celli];
				double C_SOM_temp = C_SOM[celli];

				C_DOM[celli] = kD_SOMDOM.value()/(1+kD_SOMDOM.value())*(C_DOM_temp + C_SOM_temp);
				C_SOM[celli] = 1/(1+kD_SOMDOM.value())*(C_DOM_temp + C_SOM_temp);

				C_DOM_temp = NULL;
				C_SOM_temp = NULL;

				double time_step = runTime.deltaT().value();
				celli = cellOwner[facei];
				
				//set concentration for BRNS module
				//for (int i = 0 ; i < numComp; i++){
				concentrationTemp_[0] = C_glu[celli];
				concentrationTemp_[1] = C_DOM[celli];
				concentrationTemp_[2] = 0;				//C_O2
				concentrationTemp_[3] = 0;				//C_CO2
				concentrationTemp_[4] = C_mAE[celli];
				concentrationTemp_[5] = C_mF1[celli];
				concentrationTemp_[6] = C_mF2[celli];

				boundary_flag[0] = 0;
				boundary_flag[1] = 0;
				boundary_flag[2] = 0;
				boundary_flag[3] = 0;
				boundary_flag[4] = 0;
				boundary_flag[5] = 0;
				boundary_flag[6] = 0;

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

				/*if (concentrationTemp_[2] > 1e-20)
				{ 
					C[celli] += concentrationTemp_[2];
				}
				else if ( C[celli] > 1e-20)
				{
					C[celli] += 0;
				}
				else
				{
					C[celli] = 0;
				}*/

				C_glu[celli]    = concentrationTemp_[0];//mol/m3
				C_DOM[celli]    = concentrationTemp_[1];//mol/m3
				C_O2[celli]     = concentrationTemp_[2];//mol/m3
				C_CO2[celli]    = concentrationTemp_[3];//mol/m3
				C_mAE[celli]    = concentrationTemp_[4];//mol/m3
				C_mF1[celli]    = concentrationTemp_[5];//mol/m3
				C_mF2[celli]    = concentrationTemp_[5];//mol/m3

			}
		}
	}		
}

