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
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser Base Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser Base Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Base Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef micromorphicelement_h
#define micromorphicelement_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/nlstructuralelement.h"

namespace oofem {
/**
 * Base class for second gradient continua.
 * This class contains constrained microstrain, microstretch, Couple-stress, or microdilatational continua
 * Regarding the specific formulation, new micro degrees of freedom are introduced
 * Moreover, lagrangian multiplier is introduce to constrain the micromorphic dofs to obtain the second gradient continue
 * @author Martin Horak A simplified second gradient model for dilatant materials: Theory and numerical implementation
 */
class BaseSecondGradientElement
{
protected:
    

public:
    BaseSecondGradientElement();
    virtual ~BaseSecondGradientElement() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

protected:
  
    /// Pure virtual functions
    virtual NLStructuralElement *giveElement() = 0;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &B);
    virtual void computeGmatrixAt(GaussPoint *gp, FloatMatrix &G) = 0;
  
    virtual void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeGeneralizedStressVectors(FloatArray &vP, FloatArray &vM, GaussPoint *gp, TimeStep *tStep);
    void computeDisplacementGradients(FloatArray &dudx, FloatArray &d2udx2, GaussPoint *gp, TimeStep *tStep);

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    virtual void updateInternalState(TimeStep *tStep);




};
} // end namespace oofem

#endif
