
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

IOdictionary transportProperties
(
    IOobject
    (
	"transportProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info << "Reading component C_glu" << endl;
volScalarField C_glu
(
    IOobject
    (
	"C_glu",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C_DOM" << endl;
volScalarField C_DOM
(
    IOobject
    (
	"C_DOM",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C_O2" << endl;
volScalarField C_O2
(
    IOobject
    (
	"C_O2",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C_CO2" << endl;
volScalarField C_CO2
(
    IOobject
    (
	"C_CO2",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C_SOM" << endl;
volScalarField C_SOM
(
    IOobject
    (
	"C_SOM",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C_mAE" << endl;
volScalarField C_mAE
(
    IOobject
    (
	"C_mAE",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C_mF1" << endl;
volScalarField C_mF1
(
    IOobject
    (
	"C_mF1",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading component C_mF2" << endl;
volScalarField C_mF2
(
    IOobject
    (
	"C_mF2",
	runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


dimensionedTensor D_glu
(
    transportProperties.lookup("D_glu")
);

dimensionedTensor D_DOM
(
    transportProperties.lookup("D_DOM")
);

dimensionedTensor D_O2
(
    transportProperties.lookup("D_O2")
);

dimensionedTensor D_CO2
(
    transportProperties.lookup("D_Co2")
);

dimensionedScalar kD_SOMDOM
(
    transportProperties.lookup("kD_SOMDOM")
);

/*dimensionedTensor D_mAE
(
    transportProperties.lookup("D_mAE")
);

dimensionedTensor D_mF1
(
    transportProperties.lookup("D_mF1")
);

dimensionedTensor D_mF2
(
    transportProperties.lookup("D_mF2")
);*/
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	Info << "Reading surface masters dictionary" << endl;
	wordList surfaceMastersNames;
	if (transportProperties.found("surfaceMasters"))
	{
		const dictionary& surfacesDict = transportProperties.subDict("surfaceMasters");
		Info << "Reading surface masters \n" << endl;
		surfaceMastersNames = surfacesDict.toc();
	}


    speciesTable surfaceMasters
    (
    	surfaceMastersNames
    );
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


	Info<< "Reading field Surf\n" << endl;
    volScalarField Surf
    (
        IOobject
        (
            "Surf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	//std::ostringstream oss;
	
/*	double* surfConcentration_1;
	double* surfConcentration_2;*/
	
	int ncells = mesh.cells().size();
	double* surfConcentration_2 = new double [ncells];
	double* surfArea = new double [ncells];

//	int nsurf = 2;
	
	surfArea           = (double *)malloc((size_t)(ncells * sizeof(double)));
//	surfConcentration_1  = (double *)malloc((size_t)(ncells * sizeof(double)));
//	surfConcentration_2  = (double *)malloc((size_t)(ncells * sizeof(double)));
	
	forAll(mesh.cells(),celli)
	{
		forAll(surfaceMasters, j)
		{
			double area = 0.0;
			//double mole = 0.0;
	        const volScalarField::GeometricBoundaryField& Surfbf = Surf.boundaryField();
			forAll(Surfbf,patchi)
			{
				if (Surfbf[patchi].type() == "reactingWall")
				{
					const reactingWallFvPatchScalarField& Surfcap = refCast<const reactingWallFvPatchScalarField>(Surfbf[patchi]);
					const labelList& cellOwner = Surfcap.patch().faceCells();
					const surfaceScalarField& magSf = mesh.magSf();
					const wordList& masters = Surfcap.get_surface_masters();
					const scalarList& density   = Surfcap.get_density();//mol/m^2
					forAll(masters,i)
					{
						if (masters[i]==surfaceMasters[j])
						{
							forAll(Surfbf[patchi],facei)
							{
								if (cellOwner[facei]==celli)
								{
									//mole+=density[i]*magSf.boundaryField()[patchi][facei] / mesh.V()[cellOwner[facei]];//mol/L
									area+=magSf.boundaryField()[patchi][facei] / mesh.V()[cellOwner[facei]];//m^2/m^3
								}
							}
						}
					}
				}
			}
			
			surfArea[celli]	= area;
			//Info << "surface area of cell " << celli <<"is: " << surfArea[celli] << endl;
		}
	}
	
