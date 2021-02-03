/*---------------------------------------------------------------------------*\

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

\*---------------------------------------------------------------------------*/

#include "reactingWallFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
	surfaceMasters_(dict.subDict("surfaceMasters").toc()),
	density_(surfaceMasters_.size())
{	
	forAll(surfaceMasters_, i)
	{
		const dictionary& surfaceMastersDict = dict.subDict("surfaceMasters");
		const dictionary& subdict = surfaceMastersDict.subDict(surfaceMasters_[i]);
		density_[i] = readScalar(subdict.lookup("density"));
	}
}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const reactingWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
	surfaceMasters_(ptf.surfaceMasters_),
	density_(ptf.density_)
{
}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const reactingWallFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
	surfaceMasters_(ptf.surfaceMasters_),
	density_(ptf.density_)
{}


reactingWallFvPatchScalarField::
reactingWallFvPatchScalarField
(
    const reactingWallFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
	surfaceMasters_(ptf.surfaceMasters_),
	density_(ptf.density_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reactingWallFvPatchScalarField::write(Ostream& os) const
{
	dictionary dict;

	forAll(surfaceMasters_, i)
	{
		dictionary subdict;
		subdict.add("density",density_[i]);
		dict.add(surfaceMasters_[i],subdict);
	}

	os.writeKeyword("surfaceMasters") << dict << endl;
    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    reactingWallFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
