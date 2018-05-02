/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef secondgradientms_h
#define secondgradientms_h

#include "../sm/Materials/structuralms.h"

namespace oofem {
class GaussPoint;
class Dictionary;
class Domain;

class SecondGradientMaterialStatus : public StructuralMaterialStatus
{
protected:
 
    /// Equilibrated micromorphic variable in reduced form
    FloatArray micromorphicVar;
    /// Equilibrated gradient of micromorphic variable
    FloatArray micromorphicVarGrad;
    /// Equilibrated micromorphic stress coupled to gradient of micromorphic variable
    FloatArray micromorphicStressGrad;


    /// Equilibrated micromorphic variable in reduced form
    FloatArray tempMicromorphicVar;
    /// Equilibrated gradient of micromorphic variable
    FloatArray tempMicromorphicVarGrad;
    /// Equilibrated micromorphic stress coupled to gradient of micromorphic variable
    FloatArray tempMicromorphicStressGrad;

 public:
    /// Constructor. Creates new MicromorphicMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
    SecondGradientMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~SecondGradientMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    
    /// Returns the const pointer to receiver's strain vector.
    const FloatArray &giveMicromorphicVar() const { return micromorphicVar; }
    /// Returns the const pointer to receiver's stress vector.
    const FloatArray &giveMicromorphicVarGrad() const { return micromorphicVarGrad; }
    /// Returns the const pointer to receiver's Cauchy stress vector.
    const FloatArray &giveMicromorphicStressGrad() const { return micromorphicStressGrad; }



    /// Returns the const pointer to receiver's strain vector.
    const FloatArray &giveTempMicromorphicVar() const { return tempMicromorphicVar; }
    /// Returns the const pointer to receiver's stress vector.
    const FloatArray &giveTempMicromorphicVarGrad() const { return tempMicromorphicVarGrad; }
    /// Returns the const pointer to receiver's Cauchy stress vector.
    const FloatArray &giveTempMicromorphicStressGrad() const { return tempMicromorphicStressGrad; }





    /// Assigns micromorphic variable vector to given vector mV.
    void letMicromorphicVarBe(FloatArray mv) { micromorphicVar = std :: move(mv); }
    void letMicromorphicVarGradBe(FloatArray mvg) { micromorphicVarGrad = std :: move(mvg); }

    void letMicromorphicStressGradBe(FloatArray msg) { micromorphicStressGrad = std :: move(msg); }
    /// Assigns micromorphic variable vector to given vector mV.
    void letTempMicromorphicVarBe(FloatArray mv) { tempMicromorphicVar = std :: move(mv); }
    void letTempMicromorphicVarGradBe(FloatArray mvg) { tempMicromorphicVarGrad = std :: move(mvg); }
    void letTempMicromorphicStressGradBe(FloatArray msg) { tempMicromorphicStressGrad = std :: move(msg); }
    
};
} // end namespace oofem
#endif // secondgradientms_h
