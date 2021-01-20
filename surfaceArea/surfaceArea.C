
#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #   include "setRootCase.H"

    #   include "createTime.H"
    #   include "createMesh.H"

    simpleControl simple(mesh);
	
	std::ostringstream oss;
	
	forAll(mesh.cells(),celli)
	{
		//forAll(surfaceMasters, j)
		//{
			double area = 0.0;
			//double mole = 0.0;
	        const volScalarField::GeometricBoundaryField& Surfbf = U.boundaryField();
			forAll(Surfbf,patchi)
			{
				if (Surfbf[patchi].type() == "reactingWall")
				{
					const reactingWallFvPatchScalarField& Surfcap = refCast<const reactingWallFvPatchScalarField>(Surfbf[patchi]);
					const labelList& cellOwner = Surfcap.patch().faceCells();
					const surfaceScalarField& magSf = mesh.magSf();
					//const wordList& masters = Surfcap.get_surface_masters();
					//const scalarList& density   = Surfcap.get_density();//mol/m^2
					//forAll(masters,i)
					//{
						if (masters[i]==surfaceMasters[j])
						{
							forAll(Surfbf[patchi],facei)
							{
								if (cellOwner[facei]==celli)
								{
									//mole+=density[i]*magSf.boundaryField()[patchi][facei] / mesh.V()[cellOwner[facei]];//mol/L
									area+=magSf.boundaryField()[patchi][facei]; // mesh.V()[celli][cellOwner[facei]];//m^2/m^3
								}
							}
						}
					//}
				}
			}
			
			surfArea[celli]	= area;
			//Info << "surface area of cell " << celli <<"is: " << surfArea[celli] << endl;
		//}
	}
	
	//Info <<  oss.str().c_str() << endl;
	Info <<  "area is" << area << "." << endl;
    Info<< "End\n" << endl;

    return 0;
}
