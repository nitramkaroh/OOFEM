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


#ifndef arclength_h
#define arclength_h

#include "nrsolver.h"

///@name Input fields for NRSolver
//@{
#define _IFT_ArcLength_Name "arclength"
#define  _IFT_ArcLength_numberOfAlBc "nalbc"
//@}

namespace oofem {
class Domain;
class EngngModel;

/**
 */
class OOFEM_EXPORT ArcLength : public NRSolver
{
protected:
  int numberAlBc;
  
public:
    ArcLength(Domain * d, EngngModel * m);
    virtual ~ArcLength();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "ArcLength"; }
    virtual const char *giveInputRecordName() const { return _IFT_ArcLength_Name; }

protected:
    void computeTotalLoad(FloatArray &answer, const FloatArray &R, const FloatArray &R0);
};
} // end namespace oofem
#endif // nrsolver_h
