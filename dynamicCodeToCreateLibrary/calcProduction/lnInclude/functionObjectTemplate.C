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

defineTypeNameAndDebug(calcProductionFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    calcProductionFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 0158b4e6ecd11ee2365edde8b0df226483c74425
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void calcProduction_0158b4e6ecd11ee2365edde8b0df226483c74425(bool load)
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
Foam::calcProductionFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
calcProductionFunctionObject::
calcProductionFunctionObject
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
calcProductionFunctionObject::
~calcProductionFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
calcProductionFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read calcProduction");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcProductionFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute calcProduction");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcProductionFunctionObject::write()
{
    if (false)
    {
        printMessage("write calcProduction");
    }

//{{{ begin code
    #line 28 "/home/au540668/OpenFOAM/au540668-v2412/run/WRLES/channelCases/Ret_180/casesUBulk/L2/system/controlDict/functions/calcProduction"
#include "fvCFD.H"

        // Input-Output namelist
        const word uMeanNameIn = this->dict().lookupOrDefault<word>("UMeanName", "UMean");
        const word RMeanNameIn = this->dict().lookupOrDefault<word>("RMeanName", "UPrime2Mean");
        const word prodOutNm   = this->dict().lookupOrDefault<word>("OutputName", "P_ij_RSTE");

        Info<< "Calculating Production Term for FO: " << this->name() << " at t=" << mesh().time().timeName() << endl;

        // ----- Reading Input variables ----- //

        // Lookup UMean
        if (!mesh().foundObject<volVectorField>(uMeanNameIn)) { FatalErrorInFunction << "Field " << uMeanNameIn << " not found." << exit(FatalError); }
        const volVectorField& UMean = mesh().lookupObject<volVectorField>(uMeanNameIn);
        
        // Lookup UPrime2Mean
        if (!mesh().foundObject<volSymmTensorField>(RMeanNameIn)) { FatalErrorInFunction << "Field " << RMeanNameIn << " not found." << exit(FatalError); }
        const volSymmTensorField& RMean = mesh().lookupObject<volSymmTensorField>(RMeanNameIn);
        
        Info<< "  Inputs: " << UMean.name() << ", " << RMean.name() << endl;

        // ----- Preparing Output ----- //

        // Calculating velocity gradient tensor
        volTensorField gradUMean = fvc::grad(UMean);

        // Constructing production as a volSymmTensor 
        volSymmTensorField P_ij_RSTE
        (
            IOobject(prodOutNm, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
            -2.0 * symm(RMean & gradUMean)
        );


        // Writing Output
        P_ij_RSTE.correctBoundaryConditions();
        P_ij_RSTE.write();

        Info<< "  Written field: " << P_ij_RSTE.name() << endl;
        Info<< endl;
//}}} end code

    return true;
}


bool
Foam::
calcProductionFunctionObject::end()
{
    if (false)
    {
        printMessage("end calcProduction");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

