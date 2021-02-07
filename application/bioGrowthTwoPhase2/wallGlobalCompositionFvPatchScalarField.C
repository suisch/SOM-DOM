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

#include "wallGlobalCompositionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallGlobalCompositionFvPatchScalarField::
wallGlobalCompositionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
	speciesName_(iF.name()),
	H_(0.0)
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::wallGlobalCompositionFvPatchScalarField::
wallGlobalCompositionFvPatchScalarField
(
    const wallGlobalCompositionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    speciesName_(ptf.speciesName_),
	H_(ptf.H_)
{}


Foam::wallGlobalCompositionFvPatchScalarField::
wallGlobalCompositionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
	speciesName_(iF.name()),
	H_(0.0)
{
    this->refValue() = pTraits<scalar>::zero;
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(this->refValue());
    }

    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;

	IOdictionary transportProperties 
	(
		IOobject
		(
			"transportProperties",
			this->db().time().constant(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	const dictionary& solutionSpeciesDict = transportProperties.subDict("solutionSpecies");

	if (solutionSpeciesDict.isDict(speciesName_))
    {
		const dictionary& subdict = solutionSpeciesDict.subDict(speciesName_);

    	H_ = dimensionedScalar(subdict.lookup("H")).value();
	}
}


Foam::wallGlobalCompositionFvPatchScalarField::
wallGlobalCompositionFvPatchScalarField
(
    const wallGlobalCompositionFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    speciesName_(tppsf.speciesName_),
    H_(tppsf.H_)
{}


Foam::wallGlobalCompositionFvPatchScalarField::
wallGlobalCompositionFvPatchScalarField
(
    const wallGlobalCompositionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    speciesName_(tppsf.speciesName_),
    H_(tppsf.H_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallGlobalCompositionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
	const volScalarField& alpha1 = this->db().objectRegistry::lookupObject<volScalarField>("alpha1");

	const fvPatchField<scalar>& ap = patch().patchField<volScalarField, scalar>(alpha1);
	scalarField api=ap.patchInternalField(); 

	scalarField lambda = -(1-H_)/(ap+H_*(1-ap))*(ap-api);

    this->valueFraction() = lambda/(lambda + 1);

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::wallGlobalCompositionFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wallGlobalCompositionFvPatchScalarField
    );
}

// ************************************************************************* //
