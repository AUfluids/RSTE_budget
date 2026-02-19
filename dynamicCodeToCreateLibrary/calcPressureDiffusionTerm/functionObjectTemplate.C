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

defineTypeNameAndDebug(calcPressureDiffusionTermFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    calcPressureDiffusionTermFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 85c7d85c22556eb45ea58acc8d79ed80768f4c08
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void calcPressureDiffusionTerm_85c7d85c22556eb45ea58acc8d79ed80768f4c08(bool load)
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
Foam::calcPressureDiffusionTermFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
calcPressureDiffusionTermFunctionObject::
calcPressureDiffusionTermFunctionObject
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
calcPressureDiffusionTermFunctionObject::
~calcPressureDiffusionTermFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
calcPressureDiffusionTermFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read calcPressureDiffusionTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcPressureDiffusionTermFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute calcPressureDiffusionTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcPressureDiffusionTermFunctionObject::write()
{
    if (false)
    {
        printMessage("write calcPressureDiffusionTerm");
    }

//{{{ begin code
    #line 18 "/home/au540668/OpenFOAM/au540668-v2412/run/WRLES/channelCases/Ret_180/casesUBulk/L2/system/controlDict/functions/calcPressureDiffusionTerm"
#include "fvCFD.H"

        // Field names
        const word uMeanName  = "UMean";
        const word pMeanName  = "pMean";
        const word upMeanName = "UpMean";
        const word outputDpij  = "Dp_ij_RSTE";
       
        Info<< "Calculating Pressure Diffusion Term for FO: " << this->name() << " at t=" << mesh().time().timeName() << endl;

        // Lookup Optional <u'p'>

        const volVectorField& UMean = mesh().lookupObject<volVectorField>(uMeanName);
        const volScalarField& pMean = mesh().lookupObject<volScalarField>(pMeanName);
        const volVectorField& UpMean = mesh().lookupObject<volVectorField>(upMeanName);

        volVectorField UPrimePPrimeMean (IOobject("UPrimePPrimeMean", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                         UpMean - (UMean * pMean));
        
        volTensorField gradUPrimePPrimeMean = fvc::grad(UPrimePPrimeMean); 
        
        volSymmTensorField Dp_ij_RSTE (IOobject(outputDpij, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                                       -2.0 * symm(gradUPrimePPrimeMean));
        
        Dp_ij_RSTE.correctBoundaryConditions();
        Dp_ij_RSTE.write();
            
        Info<< "  Written field: " << Dp_ij_RSTE.name() << endl;
        Info<< endl;
//}}} end code

    return true;
}


bool
Foam::
calcPressureDiffusionTermFunctionObject::end()
{
    if (false)
    {
        printMessage("end calcPressureDiffusionTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

