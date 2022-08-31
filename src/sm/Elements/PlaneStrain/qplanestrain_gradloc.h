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

#ifndef qplanestrain_gradloc_h
#define qplanestrain_gradloc_h

#include "Elements/PlaneStrain/qplanestrain.h"
#include "../sm/Elements/Micromorphic/basemicromorphicelement.h"

#define _IFT_QPlaneStrainGradLoc_Name "qplanestraingradloc"

namespace oofem {


/**
 * This class implements a quadrilateral 8-node  plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
  class QPlaneStrainGradLoc : public QPlaneStrain, public BaseMicromorphicElement
{
public:
    QPlaneStrainGradLoc(int n, Domain * d);
    virtual ~QPlaneStrainGradLoc() { }


    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QPlaneStrainGradLoc_Name; }
    virtual const char *giveClassName() const { return "QPlaneStrainGradLoc"; }

    virtual void computeMicromorphicNMatrixAt(GaussPoint *gp, FloatMatrix &N);
    virtual void computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &B);
    virtual NLStructuralElement *giveElement() { return this; }
 public:
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_m(IntArray &answer);

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep){BaseMicromorphicElement :: computeStiffnessMatrix(answer, mode, tStep);}
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord){BaseMicromorphicElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);}

    
    virtual int giveNumberOfMicromorphicDofs(){return 16;}
    virtual int giveNumberOfDisplacementDofs(){return 16;}
    virtual int giveNumberOfDofs(){return 32;}

    virtual void postInitialize();
 protected:
    virtual void computeMicromorphicVars(FloatArray &micromorphVar, FloatArray &micromorphVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep);   



};
} // end namespace oofem
#endif // qtrplanestrain_h
