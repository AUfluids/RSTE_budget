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

defineTypeNameAndDebug(calcTurbulentTransportTermFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    calcTurbulentTransportTermFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = f93081a9b9cf6aa07e2570837bb145bcfafaf89c
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void calcTurbulentTransportTerm_f93081a9b9cf6aa07e2570837bb145bcfafaf89c(bool load)
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
Foam::calcTurbulentTransportTermFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
calcTurbulentTransportTermFunctionObject::
calcTurbulentTransportTermFunctionObject
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
calcTurbulentTransportTermFunctionObject::
~calcTurbulentTransportTermFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
calcTurbulentTransportTermFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read calcTurbulentTransportTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcTurbulentTransportTermFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute calcTurbulentTransportTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcTurbulentTransportTermFunctionObject::write()
{
    if (false)
    {
        printMessage("write calcTurbulentTransportTerm");
    }

//{{{ begin code
    #line 28 "/home/au540668/OpenFOAM/au540668-v2412/run/WRLES/channelCases/Ret_180/casesUBulk/L2/system/controlDict/functions/calcTurbulentTransportTerm"
#include "fvCFD.H"

        // Field names
        const word uMeanName   = "UMean";
        const word rMeanName   = "UPrime2Mean";
        const word uuuName     = "UxUxUxMean";
        const word uuvName     = "UxUxUyMean";
        const word uuwName     = "UxUxUzMean";
        const word uvvName     = "UxUyUyMean";
        const word uvwName     = "UxUyUzMean";
        const word uwwName     = "UxUzUzMean";
        const word vvvName     = "UyUyUyMean";
        const word vvwName     = "UyUyUzMean";
        const word vwwName     = "UyUzUzMean";
        const word wwwName     = "UzUzUzMean";
        const word outputTxx   = "Txx";
        const word outputTxy   = "Txy";
        const word outputTxz   = "Txz";
        const word outputTyy   = "Tyy";
        const word outputTyz   = "Tyz";
        const word outputTzz   = "Tzz";
        const word outputTij   = "T_ij_RSTE";

        Info << "Calculating turbulent transport tensor field: " << outputTij
            << " for FO: " << this->name() << " at t=" << mesh().time().timeName() << endl;

        // Lookup input fields
        const volVectorField& UMean           = mesh().lookupObject<volVectorField>(uMeanName);
        const volSymmTensorField& UPrime2Mean = mesh().lookupObject<volSymmTensorField>(rMeanName);
        const volScalarField& UxUxUxMean      = mesh().lookupObject<volScalarField>(uuuName);
        const volScalarField& UxUxUyMean      = mesh().lookupObject<volScalarField>(uuvName);
        const volScalarField& UxUxUzMean      = mesh().lookupObject<volScalarField>(uuwName);
        const volScalarField& UxUyUyMean      = mesh().lookupObject<volScalarField>(uvvName);
        const volScalarField& UxUyUzMean      = mesh().lookupObject<volScalarField>(uvwName);
        const volScalarField& UxUzUzMean      = mesh().lookupObject<volScalarField>(uwwName);
        const volScalarField& UyUyUyMean      = mesh().lookupObject<volScalarField>(vvvName);
        const volScalarField& UyUyUzMean      = mesh().lookupObject<volScalarField>(vvwName);
        const volScalarField& UyUzUzMean      = mesh().lookupObject<volScalarField>(vwwName);
        const volScalarField& UzUzUzMean      = mesh().lookupObject<volScalarField>(wwwName);

        // Extracting components of input data
        const scalarField u   = UMean.component(vector::X);
        const scalarField v   = UMean.component(vector::Y);
        const scalarField w   = UMean.component(vector::Z);
        const scalarField uu  = UPrime2Mean.component(symmTensor::XX);
        const scalarField uv  = UPrime2Mean.component(symmTensor::XY);
        const scalarField uw  = UPrime2Mean.component(symmTensor::XZ);
        const scalarField vv  = UPrime2Mean.component(symmTensor::YY);
        const scalarField vw  = UPrime2Mean.component(symmTensor::YZ);
        const scalarField ww  = UPrime2Mean.component(symmTensor::ZZ);            
        const scalarField uuu = UxUxUxMean.internalField();
        const scalarField uuv = UxUxUyMean.internalField();
        const scalarField uuw = UxUxUzMean.internalField();
        const scalarField uvv = UxUyUyMean.internalField();
        const scalarField uvw = UxUyUzMean.internalField();
        const scalarField uww = UxUzUzMean.internalField();
        const scalarField vvv = UyUyUyMean.internalField();
        const scalarField vvw = UyUyUzMean.internalField();
        const scalarField vww = UyUzUzMean.internalField();
        const scalarField www = UzUzUzMean.internalField();
        
        // ---------- Calculate Txx ---------- //

        // Temporary vector field for calculating Txx_RSTE components
        volVectorField UxUxUk (IOobject("UxUxUk_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));
        
        // Accessing internal field
                                            
        // Calculate (<u'u'u'>, <u'u'v'>, <u'u'w'>) vector
        vectorField& resultxx = UxUxUk.ref();
        forAll(resultxx, i)
        {
            resultxx[i].x() = uuu[i] - (u[i] * uu[i] + u[i] * uu[i] + u[i] * uu[i]) - u[i]*u[i]*u[i];
            resultxx[i].y() = uuv[i] - (u[i] * uv[i] + u[i] * uv[i] + v[i] * uu[i]) - u[i]*u[i]*v[i];
            resultxx[i].z() = uuw[i] - (u[i] * uw[i] + u[i] * uw[i] + w[i] * uu[i]) - u[i]*u[i]*w[i];
        }
        UxUxUk.correctBoundaryConditions();

        // Calculate Txx by taking the divergence
        volScalarField Txx (IOobject(outputTxx, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE), 
                            fvc::div(UxUxUk));
        Txx.correctBoundaryConditions();
    
        // ---------- Calculate Txy ---------- //

        // Temporary vector field for calculating Txy_RSTE components
        volVectorField UxUyUk (IOobject("UxUyUk_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate (<u'v'u'>, <u'v'v'>, <u'v'w'>) vector    
        vectorField& resultxy = UxUyUk.ref();                                  
        forAll(resultxy, i)
        {
            resultxy[i].x() = uuv[i] - (u[i] * uv[i] + v[i] * uu[i] + u[i] * uv[i]) - u[i]*v[i]*u[i];
            resultxy[i].y() = uvv[i] - (u[i] * vv[i] + v[i] * uv[i] + v[i] * uv[i]) - u[i]*v[i]*v[i];
            resultxy[i].z() = uvw[i] - (u[i] * vw[i] + v[i] * uw[i] + w[i] * uv[i]) - u[i]*v[i]*w[i];
        }
        UxUyUk.correctBoundaryConditions();

        // Calculate Txy by taking the divergence            
        volScalarField Txy (IOobject(outputTxy, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(UxUyUk));
        Txy.correctBoundaryConditions();
    
        // ---------- Calculate Txz ---------- //

        // Temporary vector field for calculating Txz_RSTE components
        volVectorField UxUzUk (IOobject("UxUzUk_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate (<u'w'u'>, <u'w'v'>, <u'w'w'>) vector
        vectorField& resultxz = UxUzUk.ref();
        forAll(resultxz, i)
        {
            resultxz[i].x() = uuw[i] - (u[i] * uw[i] + w[i] * uu[i] + u[i] * uw[i]) - u[i]*w[i]*u[i];
            resultxz[i].y() = uvv[i] - (u[i] * vw[i] + w[i] * uv[i] + v[i] * uw[i]) - u[i]*w[i]*v[i];
            resultxz[i].z() = uww[i] - (u[i] * ww[i] + w[i] * uw[i] + w[i] * uw[i]) - u[i]*w[i]*w[i];
        }
        UxUzUk.correctBoundaryConditions();

        // Calculate Txz by taking the divergence            
        volScalarField Txz (IOobject(outputTxz, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(UxUzUk));
        Txz.correctBoundaryConditions();
    
        // ---------- Calculate Tyy ---------- //  
            
        // Temporary vector field for calculating Tyy_RSTE components
        volVectorField UyUyUk (IOobject("UyUyUk_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate (<v'v'u'>, <v'v'v'>, <v'v'w'>) vector                                  
        vectorField& resultyy = UyUyUk.ref();                                  
        forAll(resultyy, i)
        {
            resultyy[i].x() = uvv[i] - (v[i] * uv[i] + v[i] * uv[i] + u[i] * vv[i]) - v[i]*v[i]*u[i];
            resultyy[i].y() = vvv[i] - (v[i] * vv[i] + v[i] * vv[i] + v[i] * vv[i]) - v[i]*v[i]*v[i];
            resultyy[i].z() = vvw[i] - (v[i] * vw[i] + v[i] * vw[i] + w[i] * vv[i]) - v[i]*v[i]*w[i];
        }
        UyUyUk.correctBoundaryConditions();

        // Calculate Tyy by taking the divergence                        
        volScalarField Tyy (IOobject(outputTyy, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(UyUyUk));
        Tyy.correctBoundaryConditions();
        
        // ---------- Calculate Tyz ---------- //  
            
        // Temporary vector field for calculating Tyz_RSTE components
        volVectorField UyUzUk (IOobject("UyUzUk_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));
        

        // Calculate (<v'w'u'>, <v'w'v'>, <v'w'w'>) vector 
        vectorField& resultyz = UyUzUk.ref();            
        forAll(resultyz, i)
        {
            resultyz[i].x() = uvw[i] - (v[i] * uw[i] + w[i] * uv[i] + u[i] * vw[i]) - v[i]*w[i]*u[i];
            resultyz[i].y() = vvw[i] - (v[i] * vw[i] + w[i] * vv[i] + v[i] * vw[i]) - v[i]*w[i]*v[i];
            resultyz[i].z() = vww[i] - (v[i] * ww[i] + w[i] * vw[i] + w[i] * vw[i]) - v[i]*w[i]*w[i];
        }
        UyUzUk.correctBoundaryConditions();

        // Calculate Tyz by taking the divergence                                    
        volScalarField Tyz (IOobject(outputTyz, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(UyUzUk));
        Tyz.correctBoundaryConditions();

        // ---------- Calculate Tzz ---------- //               

        // Temporary vector field for calculating Tzz_RSTE components
        volVectorField UzUzUk (IOobject("UzUzUk_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));


        // Calculate (<w'w'u'>, <w'w'v'>, <w'w'w'>) vector             
        vectorField& resultzz = UzUzUk.ref();
        forAll(resultzz, i)
        {
            resultzz[i].x() = uww[i] - (w[i] * uw[i] + w[i] * uw[i] + u[i] * ww[i]) - w[i]*w[i]*u[i];
            resultzz[i].y() = vww[i] - (w[i] * vw[i] + w[i] * vw[i] + v[i] * ww[i]) - w[i]*w[i]*v[i];
            resultzz[i].z() = www[i] - (w[i] * ww[i] + w[i] * ww[i] + w[i] * ww[i]) - w[i]*w[i]*w[i];
        }
        UzUzUk.correctBoundaryConditions();

        // Calculate Tzz by taking the divergence                                                
        volScalarField Tzz (IOobject(outputTzz, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(UzUzUk));
        Tzz.correctBoundaryConditions();

        // ---------- Build full symmTensor ---------- //  
        
        volSymmTensorField Tij (IOobject(outputTij, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                                mesh(),dimensionedSymmTensor("zero", Txx.dimensions(), symmTensor::zero));

        symmTensorField& TijField = Tij.ref();

        forAll(TijField, i)
        {
            TijField[i].xx() = -Txx[i];
            TijField[i].xy() = -Txy[i];
            TijField[i].xz() = -Txz[i];
            TijField[i].yy() = -Tyy[i];
            TijField[i].yz() = -Tyz[i];
            TijField[i].zz() = -Tzz[i];
        }

        Tij.correctBoundaryConditions();
        Tij.write();

        Info << "  Written symmTensor field: " << Tij.name() << endl;
        Info<< endl;
//}}} end code

    return true;
}


bool
Foam::
calcTurbulentTransportTermFunctionObject::end()
{
    if (false)
    {
        printMessage("end calcTurbulentTransportTerm");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

