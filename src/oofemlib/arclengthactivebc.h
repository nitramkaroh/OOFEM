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

#ifndef arclengthactivebc_h
#define arclengthactivebc_h

#include "activebc.h"
#include "floatmatrix.h"


namespace oofem {
/**
 * Abstract base class for arc lenght active boundary conditions.
 * This boundary condition specify the arc length constrain, i.e., it can have the classical form
 * (\Delta u)^T\Delta u  + \Delta \lambda \Psi \bar{f}_{ext}^t\bar{f}_{ext} = l_{arc}
 * or another form, like the dissipation based one, i.e., \Delta D_{diss} = l_{diss}
 */
class OOFEM_EXPORT ArcLengthActiveBoundaryCondition : public ActiveBoundaryCondition
{
 protected:
  FloatMatrix referenceLoadMatrix;
  double lambda, deltaLambda;
  
public:
    /**
     * Constructor. Creates boundary and active condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    ArcLengthActiveBoundaryCondition(int n, Domain * d) : ActiveBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~ArcLengthActiveBoundaryCondition() { }
    virtual double giveDeltaLambda() = 0;
    virtual void setReferenceLoadVector(const FloatArray &rF){referenceLoadMatrix.initFromVector(rF, false);}


};
} // end namespace oofem
#endif // arclengthactivebc_h
