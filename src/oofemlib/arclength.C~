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

#include "arclength.h"
#include "timestep.h"
#include "mathfem.h"
#include "activebc.h"
#include "function.h"
#include "linesearch.h"
#include "classfactory.h"
#include "exportmodulemanager.h"
#include "engngm.h"
#include "parallelcontext.h"
#include "unknownnumberingscheme.h"

#ifdef __PETSC_MODULE
 #include "petscsolver.h"
 #include "petscsparsemtrx.h"
#endif

#include <cstdio>




namespace oofem {

REGISTER_SparseNonLinearSystemNM(ArcLength)
  
ArcLength :: ArcLength(Domain *d, EngngModel *m) : NRSolver(d, m)
{
  
}


ArcLength :: ~ArcLength()
{
}


IRResultType
ArcLength :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    result  = NRSolver :: initializeFrom(ir);      
    IR_GIVE_FIELD(ir, numberAlBc, _IFT_ArcLength_numberOfAlBc);
    //arcLengthBC = engngm->giveBoundaryCondition(numberAlBc);
    //dynamic cast ...
    return result;
}

void
ArcLength :: computeTotalLoad(FloatArray &answer, const FloatArray &R, const FloatArray &R0)
{
  //double lambda = arcLengthBC->giveLambda();
  double lambda = 1;
  answer = R;
  answer.times(lambda);
  if ( R0.giveSize() ) {
    answer.add(R0);
  }
}


  


} // end namespace oofem
