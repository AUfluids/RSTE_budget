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

defineTypeNameAndDebug(calcPressureStrainTermFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    calcPressureStrainTermFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 4792976e08e7ed5d145394b349c344fc7714027b
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void calcPressureStrainTerm_4792976e08e7ed5d145394b349c344fc7714027b(bool load)
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
Foam::calcPressureStrainTermFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
calcPressureStrainTermFunctionObject::
calcPressureStrainTermFunctionObject
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
calcPressureStrainTermFunctionObject::
~calcPressureStrainTermFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
calcPressureStrainTermFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read calcPressureStrainTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcPressureStrainTermFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute calcPressureStrainTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcPressureStrainTermFunctionObject::write()
{
    if (false)
    {
        printMessage("write calcPressureStrainTerm");
    }

//{{{ begin code
    #line 18 "/home/au540668/OpenFOAM/au540668-v2412/run/WRLES/channelCases/Ret_180/casesUBulk/L2/system/controlDict/functions/calcPressureStrainTerm"
#include "fvCFD.H"

        // Field names
        const word UMeanName      = "UMean";
        const word pMeanName      = "pMean";
        const word pGradUMeanName = "pGradUMean";
        const word outputPhiij    = "Phi_ij_RSTE";        
        
        Info<< "Calculating Pressure-Strain (from means) Term for FO: " << this->name() << " at t=" << mesh().time().timeName() << endl;

        // Loading UMean
        if (!mesh().foundObject<volVectorField>(UMeanName))
        {
            FatalErrorInFunction << "Input field " << UMeanName << " not found for FO " << this->name() << "."
                << exit(FatalError);
        }
        const volVectorField& UMean = mesh().lookupObject<volVectorField>(UMeanName);
        Info<< "  Using input UMean field: " << UMean.name() << endl;

        // Loading pMean
        if (!mesh().foundObject<volScalarField>(pMeanName))
        {
            FatalErrorInFunction << "Input field " << pMeanName << " not found for FO " << this->name() << "."
                << exit(FatalError);
        }
        const volScalarField& pMean = mesh().lookupObject<volScalarField>(pMeanName);
        Info<< "  Using input pMean field: " << pMean.name() << endl;

        // Loading pGradUMean
        if (!mesh().foundObject<volTensorField>(pGradUMeanName))
        {
            FatalErrorInFunction << "Input field " << pGradUMeanName << " not found for FO " << this->name() << "."
                << exit(FatalError);
        }
        const volTensorField& pGradUMean = mesh().lookupObject<volTensorField>(pGradUMeanName);
        Info<< "  Using input pGradUMean field: " << pGradUMean.name() << endl;


        Info<< "  Calculating Sij_mean = symm(fvc::grad(UMean))" << endl;
        volSymmTensorField Sij_mean = symm(fvc::grad(UMean));

        Info<< "  Calculating pSij_mean = symm(pGradUMean)" << endl;
        volSymmTensorField pSij_mean = symm(pGradUMean);

        Info<< "  Calculating " << outputPhiij << " = pMean * Sij_mean" << endl;
        volSymmTensorField Phiij (IOobject(outputPhiij, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
            2.0 * (pSij_mean - pMean * Sij_mean));

        Phiij.correctBoundaryConditions();
        Phiij.write();
    
        Info << "  Written symmTensor field: " << Phiij.name() << endl;
        Info<< endl;
//}}} end code

    return true;
}


bool
Foam::
calcPressureStrainTermFunctionObject::end()
{
    if (false)
    {
        printMessage("end calcPressureStrainTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

