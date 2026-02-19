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

defineTypeNameAndDebug(calcConvectionFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    calcConvectionFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = f80c8ad2b64925e9fc32c424e676a62d73cb56b8
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void calcConvection_f80c8ad2b64925e9fc32c424e676a62d73cb56b8(bool load)
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
Foam::calcConvectionFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
calcConvectionFunctionObject::
calcConvectionFunctionObject
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
calcConvectionFunctionObject::
~calcConvectionFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
calcConvectionFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read calcConvection");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcConvectionFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute calcConvection");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcConvectionFunctionObject::write()
{
    if (false)
    {
        printMessage("write calcConvection");
    }

//{{{ begin code
    #line 72 "/home/au540668/OpenFOAM/au540668-v2412/run/WRLES/channelCases/Ret_180/casesUBulk/L2/system/controlDict/functions/calcConvection"
// Importing ...  
        #include "fvCFD.H"

        // Field names
        const word uMeanName   = "UMean";
        const word rMeanName   = "UPrime2Mean";
        const word outputTxx   = "C_xx";
        const word outputTxy   = "C_xy";
        const word outputTxz   = "C_xz";
        const word outputTyy   = "C_yy";
        const word outputTyz   = "C_yz";
        const word outputTzz   = "C_zz";
        const word output_C_ij = "C_ij_RSTE";

        Info << "Calculating convection tensor field: " << output_C_ij
            << " for FO: " << this->name() << " at t=" << mesh().time().timeName() << endl;

        // Lookup input fields
        const volVectorField& UMean           = mesh().lookupObject<volVectorField>(uMeanName);
        const volSymmTensorField& UPrime2Mean = mesh().lookupObject<volSymmTensorField>(rMeanName);


        // Extracting components of input data
        const scalarField u  = UMean.component(vector::X);
        const scalarField v  = UMean.component(vector::Y);
        const scalarField w  = UMean.component(vector::Z);
        const scalarField uu = UPrime2Mean.component(symmTensor::XX);
        const scalarField uv  = UPrime2Mean.component(symmTensor::XY);
        const scalarField uw  = UPrime2Mean.component(symmTensor::XZ);
        const scalarField vv  = UPrime2Mean.component(symmTensor::YY);
        const scalarField vw  = UPrime2Mean.component(symmTensor::YZ);
        const scalarField ww  = UPrime2Mean.component(symmTensor::ZZ);

        // ---------- Calculate xx component ---------- //               

        volVectorField u_uu (IOobject("u_uu_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate vector...             
        vectorField& resultxx = u_uu.ref();
        forAll(resultxx, i)
        {
            resultxx[i].x() = u[i] * uu[i];
            resultxx[i].y() = v[i] * uu[i];
            resultxx[i].z() = w[i] * uu[i];
        }
        u_uu.correctBoundaryConditions();
        
        // Calculate divergence...                                                
        volScalarField C_xx (IOobject("C_xx", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(u_uu));
        C_xx.correctBoundaryConditions();

        // ---------- Calculate xy component ---------- //               

        volVectorField u_uv (IOobject("u_uv_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate vector...             
        vectorField& resultxy = u_uv.ref();
        forAll(resultxy, i)
        {
            resultxy[i].x() = u[i] * uv[i];
            resultxy[i].y() = v[i] * uv[i];
            resultxy[i].z() = w[i] * uv[i];
        }
        u_uv.correctBoundaryConditions();
        
        // Calculate divergence...                                                
        volScalarField C_xy (IOobject("C_xy", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(u_uv));
        C_xy.correctBoundaryConditions();

        // ---------- Calculate xz component ---------- //               

        volVectorField u_uw (IOobject("u_uw_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate vector...             
        vectorField& resultxz = u_uw.ref();
        forAll(resultxz, i)
        {
            resultxz[i].x() = u[i] * uw[i];
            resultxz[i].y() = v[i] * uw[i];
            resultxz[i].z() = w[i] * uw[i];
        }
        u_uw.correctBoundaryConditions();
        
        // Calculate divergence...                                                
        volScalarField C_xz (IOobject("C_xz", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(u_uw));
        C_xz.correctBoundaryConditions();        

        // ---------- Calculate yy component ---------- //               

        volVectorField u_vv (IOobject("u_vv_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate vector...             
        vectorField& resultyy = u_vv.ref();
        forAll(resultyy, i)
        {
            resultyy[i].x() = u[i] * vv[i];
            resultyy[i].y() = v[i] * vv[i];
            resultyy[i].z() = w[i] * vv[i];
        }
        u_vv.correctBoundaryConditions();
        
        // Calculate divergence...                                                
        volScalarField C_yy (IOobject("C_yy", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(u_vv));
        C_yy.correctBoundaryConditions();        

        // ---------- Calculate yz component ---------- //               

        volVectorField u_vw (IOobject("u_vw_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate vector...             
        vectorField& resultyz = u_vw.ref();
        forAll(resultyz, i)
        {
            resultyz[i].x() = u[i] * vw[i];
            resultyz[i].y() = v[i] * vw[i];
            resultyz[i].z() = w[i] * vw[i];
        }
        u_vw.correctBoundaryConditions();
        
        // Calculate divergence...                                                
        volScalarField C_yz (IOobject("C_yz", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(u_vw));
        C_yz.correctBoundaryConditions();                
    
        // ---------- Calculate zz component ---------- //               

        volVectorField u_ww (IOobject("u_ww_tmp", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                                mesh(), dimensionedVector("zero", dimVelocity*dimVelocity*dimVelocity, vector::zero));

        // Calculate vector...             
        vectorField& resultzz = u_ww.ref();
        forAll(resultzz, i)
        {
            resultzz[i].x() = u[i] * ww[i];
            resultzz[i].y() = v[i] * ww[i];
            resultzz[i].z() = w[i] * ww[i];
        }
        u_ww.correctBoundaryConditions();
        
        // Calculate divergence...                                                
        volScalarField C_zz (IOobject("C_zz", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                            fvc::div(u_ww));
        C_zz.correctBoundaryConditions();                

        // ---------- Build full symmTensor ---------- //  
        
        volSymmTensorField C_ij(IOobject(output_C_ij, mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
                                   mesh(), dimensionedSymmTensor("zero", dimensionSet(0, 1, -2, 0, 0, 0, 0), symmTensor::zero));

        symmTensorField& C_ij_Field = C_ij.ref();

        forAll(C_ij_Field, i)
        {
            C_ij_Field[i].xx() = C_xx[i];
            C_ij_Field[i].xy() = C_xy[i];
            C_ij_Field[i].xz() = C_xz[i];
            C_ij_Field[i].yy() = C_yy[i];
            C_ij_Field[i].yz() = C_yz[i];
            C_ij_Field[i].zz() = C_zz[i];
        }

        C_ij.correctBoundaryConditions();
        C_ij.write();

        Info << "  Written symmTensor field: " << C_ij.name() << endl;
        Info<< endl;
//}}} end code

    return true;
}


bool
Foam::
calcConvectionFunctionObject::end()
{
    if (false)
    {
        printMessage("end calcConvection");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

