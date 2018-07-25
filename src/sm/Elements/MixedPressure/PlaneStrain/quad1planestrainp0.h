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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#ifndef quad1planestrainp0_h
#define quad1planestrainp0_h

#include "../sm/Elements/PlaneStrain/quad1planestrain.h"
#include "../sm/Elements/Micromorphic/basemixedpressureelement.h"

#define _IFT_PlaneStrain_Name "planestrainmicropolar"

namespace oofem {
class FEI2dQuadLin;

class Quad1PlaneStrainP0 : public Quad1PlaneStrain, public BaseMicromorphicElement
{
protected:
    static FEI2dQuadLin interpolation;

public:
    Quad1PlaneStrainP0(int n, Domain * d);
    virtual ~Quad1PlaneStrainP0() { }

 protected:


    virtual void computePressureNMatrixAt(GaussPoint *gp, FloatMatrix &N);
    virtual void computeDeviatoricVolumetricBMatrices(FloatMatrix &Bdev, FloatMatrix &Bvol, GaussPoint *gp, NLStructuralElement *element);

    virtual NLStructuralElement *giveElement() { return this; }
 public:
    virtual const char *giveInputRecordName() const { return _IFT_Quad1PlaneStrainP0_Name; }
    virtual const char *giveClassName() const { return "Quad1PlaneStrainP0"; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_m(IntArray &answer);

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep){BaseMixedPressureElement :: computeStiffnessMatrix(answer, mode, tStep);}
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord){BaseMixedPressureElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);}


    virtual int giveNumberOfPressureDofs(){return 4;}
    virtual int giveNumberOfDisplacementDofs(){return 8;}
    virtual int giveNumberOfDofs(){return 12;}

};
} // end namespace oofem
#endif // quad1planestrainp0_h
