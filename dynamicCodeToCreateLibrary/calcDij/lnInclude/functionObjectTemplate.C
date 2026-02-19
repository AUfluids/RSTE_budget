/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
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

#include "functionObjectTemplate.H"
#define namespaceFoam  // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(calcDijFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    calcDijFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 5485f91b04b9ad1e58bf5bbbdeb48ad89875cbec
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void calcDij_5485f91b04b9ad1e58bf5bbbdeb48ad89875cbec(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMesh&
Foam::calcDijFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
calcDijFunctionObject::
calcDijFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
calcDijFunctionObject::
~calcDijFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
calcDijFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read calcDij");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcDijFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute calcDij");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcDijFunctionObject::write()
{
    if (false)
    {
        printMessage("write calcDij");
    }

//{{{ begin code
    #line 28 "/home/au540668/OpenFOAM/au540668-v2412/run/WRLES/channelCases/Ret_180/casesUBulk/L2/system/controlDict/functions/calcViscousDiffusionTerm"
#include "fvCFD.H"

        const word RMeanNameIn      = this->dict().lookupOrDefault<word>("RMeanName", "UPrime2Mean");
        const word transportPropsNm = this->dict().lookupOrDefault<word>("transportProperties", "transportProperties");
        const word viscDiffOutNm    = this->dict().lookupOrDefault<word>("OutputName", "D_ij_RSTE");

        Info<< "Calculating Viscous Diffusion Term for FO: " << this->name() << " at t=" << mesh().time().timeName() << endl;

        // Get nu
        if (!mesh().foundObject<IOdictionary>(transportPropsNm))
        {
            FatalErrorInFunction << "Transport properties '" << transportPropsNm << "' not found." << exit(FatalError);
        }
        const dictionary& transportProps = mesh().lookupObject<IOdictionary>(transportPropsNm);
        dimensionedScalar nu("nu", dimViscosity, transportProps);

        // Lookup RMean
        if (!mesh().foundObject<volSymmTensorField>(RMeanNameIn))
        {
            FatalErrorInFunction << "Field " << RMeanNameIn << " not found." << exit(FatalError);
        }
        const volSymmTensorField& RMean = mesh().lookupObject<volSymmTensorField>(RMeanNameIn);
        Info<< "  Inputs: " << RMean.name() << ", nu=" << nu << endl;

        volSymmTensorField D_ij_RSTE (
            IOobject(viscDiffOutNm, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
            nu * fvc::laplacian(RMean)
        );

        const vectorField& cellCentres = mesh().C();
        const vectorField& faceCentres = mesh().Cf();

        forAll(mesh().boundary(), patchI)
        {
            const fvPatch& patch = mesh().boundary()[patchI];

            if (isA<wallFvPatch>(patch))
            {
                const label startFaceI = patch.start();
                const labelList& faceCells = patch.faceCells();

                forAll(faceCells, faceI)
                {
                    const label c0 = faceCells[faceI];
                    const label meshFaceI = startFaceI + faceI;

                    const point faceCenter = faceCentres[meshFaceI];
                    const vector normal = (cellCentres[c0] - faceCenter) / mag(cellCentres[c0] - faceCenter);

                    const labelList& nbrs = mesh().cellCells()[c0];
                    //Info << "Cell " << c0 << " (" << patch.name() << ") has " << nbrs.size() << " neighbors." << endl;

                    label bestNbr = -1;
                    scalar maxProj = -GREAT;

                    for (const label nbr : nbrs)
                    {
                        const vector d = cellCentres[nbr] - faceCenter;
                        scalar proj = d & normal;

                        //Info << "    Neighbor " << nbr << " projection: " << proj << endl;

                        if (proj > maxProj)
                        {
                            maxProj = proj;
                            bestNbr = nbr;
                        }
                    }

                    if (bestNbr != -1)
                    {
                        D_ij_RSTE[c0] = D_ij_RSTE[bestNbr];
                        //Info << "  Overwriting cell " << c0
                        //     << " from neighbor " << bestNbr
                        //     << " with projection " << maxProj << endl;
                    }
                    else
                    {
                        //Info << "  WARNING: No valid neighbor found for cell " << c0
                        //     << " on patch " << patch.name() << endl;
                    }
                }
            }
        }

        D_ij_RSTE.correctBoundaryConditions();
        D_ij_RSTE.write();

        Info<< "  Written field: " << D_ij_RSTE.name() << endl;
        Info<< endl;
//}}} end code

    return true;
}


bool
Foam::
calcDijFunctionObject::end()
{
    if (false)
    {
        printMessage("end calcDij");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

