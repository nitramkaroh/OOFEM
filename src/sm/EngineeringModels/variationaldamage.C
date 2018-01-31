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

#include "../sm/EngineeringModels/variationaldamage.h"
#include "../sm/Elements/structuralelement.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "error.h"
#include "verbose.h"
#include "sparsenonlinsystemnm.h"
#include "nrsolver.h"
#include "staggeredsolver.h"
#include "calmls.h"
#include "classfactory.h"
#include "sparsemtrx.h"
#include "mathfem.h"
#include "dofmanager.h"
#include "dof.h"
#include "unknownnumberingscheme.h"
#include "activebc.h"
#include "set.h"
#include "load.h"
#include "bodyload.h"
#include "boundaryload.h"
#include "nodalload.h"
#include "feinterpol3d.h"




namespace oofem {
REGISTER_EngngModel(VariationalDamage);

VariationalDamage :: VariationalDamage(int i, EngngModel *_master) : NonLinearStatic(i, _master),    damageIndicatorArray()
{

}


VariationalDamage :: ~VariationalDamage()
{
}




IRResultType
VariationalDamage :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = NonLinearStatic :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    IR_GIVE_FIELD(ir, this->totalIdList, _IFT_VariationalDamage_DofIdList);
    IR_GIVE_FIELD(ir, this->idPos, _IFT_VariationalDamage_DofIdListPositions);
    this->maxActivNodes = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->maxActivNodes, _IFT_VariationalDamage_MaxActivatedNodes);

    //    this->instanciateYourself();
    return IRRT_OK;
}


  /*NumericalMethod*
VariationalDamage :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
      //nMethod =  new StaggeredSolver(this->giveDomain(1), this);
      nMethod =  new NRSolver(this->giveDomain(1), this);
    }
    return this->nMethod;
    }*/


void
VariationalDamage :: updateYourself(TimeStep *tStep)
{
    StructuralEngngModel :: updateYourself(tStep);
    damageIndicatorArray.clear();
}

  
void
VariationalDamage :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
//
// updates some component, which is used by numerical method
// to newly reached state. used mainly by numerical method
// when new tangent stiffness is needed during finding
// of new equilibrium stage.
//
{
    switch ( cmpn ) {
    case NonLinearLhs:
      stiffnessMatrix->zero(); // zero stiffness matrix
#ifdef VERBOSE
      OOFEM_LOG_DEBUG("Assembling tangent stiffness matrix\n");
#endif
      this->assemble(*stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
		     EModelDefaultEquationNumbering(), d);       
        break;
    case InternalRhs:
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Updating internal forces\n");
#endif
        // update internalForces and internalForcesEBENorm concurrently
        this->giveInternalForces(internalForces, true, d->giveNumber(), tStep);
        break;
    case InitialGuess:
      
      this-> giveInitialGuess(d->giveNumber(), tStep);
      break;
    default:
        OOFEM_ERROR("Unknown Type of component.");
    }
}





void 
VariationalDamage :: initializeYourself(TimeStep *tStep)
{
    NonLinearStatic :: initializeYourself(tStep);
    int numDofIdGroups = idPos.giveSize()/2;
    this->UnknownNumberingSchemeList.resize(numDofIdGroups);
    IntArray idList; 
    for ( int i = 0; i < numDofIdGroups; i++ ) {
        int sz = idPos.at(i*2 + 2) - idPos.at(i*2 + 1) + 1;
        idList.resize(sz);
        for ( int j = 1; j <= sz; j++) {
            int pos = idPos.at(i*2 + 1) + j - 1;
            idList.at(j) = totalIdList.at(pos);
        }
        this->UnknownNumberingSchemeList[i].setDofIdArray(idList);
        this->UnknownNumberingSchemeList[i].setNumber(i+1);
  
    }    


    this->locArrayList.resize(numDofIdGroups);            
    
    for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
      this->giveTotalLocationArray(this->locArrayList[dG], UnknownNumberingSchemeList[dG], this->giveDomain(1));      }
}


void
VariationalDamage :: giveTotalLocationArray(IntArray &condensedLocationArray, const UnknownNumberingScheme &s, Domain *d)
{
    IntArray nodalArray, ids, locationArray;
    locationArray.clear();
    
    for ( auto &dman : d->giveDofManagers() ) {
        dman->giveCompleteLocationArray(nodalArray, s);
        locationArray.followedBy(nodalArray);
    }
    for ( auto &elem : d->giveElements() ) {
        for ( int i = 1; i <= elem->giveNumberOfInternalDofManagers(); i++ ) {
            elem->giveInternalDofManDofIDMask(i, ids);
            elem->giveInternalDofManager(i)->giveLocationArray(ids, nodalArray, s);
            locationArray.followedBy(nodalArray);
        }
    }
    
    
    IntArray nonZeroMask;
    nonZeroMask.findNonzeros(locationArray);

    condensedLocationArray.resize(nonZeroMask.giveSize());
    for ( int i = 1; i <= nonZeroMask.giveSize(); i++ ) {
        condensedLocationArray.at(i) = locationArray.at( nonZeroMask.at(i) );    
    }
}    

  
  

void
VariationalDamage :: giveInternalForces(FloatArray &answer, bool normFlag, int di, TimeStep *tStep)
{

    // Simply assembles contributions from each element in domain
    Domain *domain = this->giveDomain(di);
    // Update solution state counter
    tStep->incrementStateCounter();

    answer.resize( this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() ) );
    answer.zero();
    this->assembleVector(answer, tStep, InternalForceAssembler(), VM_Total,
                         EModelDefaultEquationNumbering(), domain, normFlag ? & this->internalForcesEBENorm : NULL);

    // Redistributes answer so that every process have the full values on all shared equations
    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);

    // Remember last internal vars update time stamp.
    internalVarUpdateStamp = tStep->giveSolutionStateCounter();

    FloatArray iFdamage;
    iFdamage.beSubArrayOf(answer, locArrayList[1] );        

    damageIndicatorArray.clear();

    FloatArray copyiFd(iFdamage);
    double min = 0;
    if(maxActivNodes > 0) {
      // find nodes which should be activated
      for(int i = 1; i <= maxActivNodes; i++) {	  
	int minInd = copyiFd.giveIndexMinElem();
	double newMin = copyiFd.at(minInd);
	if(newMin > 0) {
	  break;
	} else {
	  min = newMin;
	}

	copyiFd.at(minInd) = 1;
	iFdamage.at(minInd) *= 2;
      }
    }
      
    // find nodes which should be activated
    for(int i = 1; i <= iFdamage.giveSize(); i++) {	  
      if(iFdamage.at(i) >= min) {	 
	answer.at(locArrayList[1].at(i)) = 0;
      } else {
	damageIndicatorArray.followedBy(locArrayList[1].at(i));
      }
    }



    
}


void VariationalDamage :: assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
                            const UnknownNumberingScheme &s, Domain *domain)
//
// assembles matrix
//
{
    IntArray loc;
    FloatMatrix mat, R;

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    int nelem = domain->giveNumberOfElements();
#ifdef _OPENMP
 #pragma omp parallel for shared(answer) private(mat, R, loc)
#endif
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = domain->giveElement(ielem);
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote || !element->isActivated(tStep) ) {
            continue;
        }

        ma.matrixFromElement(mat, *element, tStep);

        if ( mat.isNotEmpty() ) {
            ma.locationFromElement(loc, *element, s);
            ///@todo This rotation matrix is not flexible enough.. it can only work with full size matrices and doesn't allow for flexibility in the matrixassembler.
            if ( element->giveRotationMatrix(R) ) {
                mat.rotatedWith(R);
            }

#ifdef _OPENMP
 #pragma omp critical
#endif

	    for(int i = 1; i<= loc.giveSize(); i++) {
	      int index = damageIndicatorArray.findFirstIndexOf(loc.at(i));
	      int damageLoc = locArrayList[1].findFirstIndexOf(loc.at(i));
	      /*	      if(index != 0) {
		FloatMatrix m = {{1}};
		IntArray ll = {loc.at(i)};
		loc.at(i) = 0;
		answer.assemble(ll, m);
		}*/
	      if(index == 0 && damageLoc != 0) {
		FloatMatrix m = {{1}};
		IntArray ll = {loc.at(i)};
		loc.at(i) = 0;
		answer.assemble(ll, m);
	      }

	      
	    }
	    

	    
            if ( answer.assemble(loc, mat) == 0 ) {
                OOFEM_ERROR("sparse matrix assemble error");
            }
        }
    }

    int nbc = domain->giveNumberOfBoundaryConditions();
    for ( int i = 1; i <= nbc; ++i ) {
        GeneralBoundaryCondition *bc = domain->giveBc(i);
        ActiveBoundaryCondition *abc;
        Load *load;

        if ( ( abc = dynamic_cast< ActiveBoundaryCondition * >(bc) ) ) {
            ma.assembleFromActiveBC(answer, *abc, tStep, s, s);
        } else if ( bc->giveSetNumber() && ( load = dynamic_cast< Load * >(bc) ) && bc->isImposed(tStep) ) {
            // Now we assemble the corresponding load type for the respective components in the set:
            IntArray loc, bNodes;
            FloatMatrix mat, R;
            BodyLoad *bodyLoad;
            SurfaceLoad* sLoad;
            EdgeLoad* eLoad;
            Set *set = domain->giveSet( bc->giveSetNumber() );

            if ( ( bodyLoad = dynamic_cast< BodyLoad * >(load) ) ) { // Body load:
                const IntArray &elements = set->giveElementList();
                for ( int ielem = 1; ielem <= elements.giveSize(); ++ielem ) {
                    Element *element = domain->giveElement( elements.at(ielem) );
                    mat.clear();
                    ma.matrixFromLoad(mat, *element, bodyLoad, tStep);

                    if ( mat.isNotEmpty() ) {
                        if ( element->giveRotationMatrix(R) ) {
                            mat.rotatedWith(R);
                        }

                        ma.locationFromElement(loc, *element, s);
                        answer.assemble(loc, mat);
                    }
                }
            } else if ( ( sLoad = dynamic_cast< SurfaceLoad * >(load) ) ) { // Surface load:
                const IntArray &boundaries = set->giveBoundaryList();
                for ( int ibnd = 1; ibnd <= boundaries.giveSize() / 2; ++ibnd ) {
                    Element *element = domain->giveElement( boundaries.at(ibnd * 2 - 1) );
                    int boundary = boundaries.at(ibnd * 2);
                    mat.clear();
                    ma.matrixFromSurfaceLoad(mat, *element, sLoad, boundary, tStep);

                    if ( mat.isNotEmpty() ) {
                        element->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
                        if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                            mat.rotatedWith(R);
                        }

                        ma.locationFromElementNodes(loc, *element, bNodes, s);
                        answer.assemble(loc, mat);
                    }
                }
            } else if ( ( eLoad = dynamic_cast< EdgeLoad * >(load) ) ) { // Edge load:
                const IntArray &edgeBoundaries = set->giveEdgeList();
                for ( int ibnd = 1; ibnd <= edgeBoundaries.giveSize() / 2; ++ibnd ) {
                    Element *element = domain->giveElement( edgeBoundaries.at(ibnd * 2 - 1) );
                    int boundary = edgeBoundaries.at(ibnd * 2);
                    mat.clear();
                    ma.matrixFromEdgeLoad(mat, *element, eLoad, boundary, tStep);

                    if ( mat.isNotEmpty() ) {
                        element->giveInterpolation()->boundaryEdgeGiveNodes(bNodes, boundary);
                        if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                            mat.rotatedWith(R);
                        }

                        ma.locationFromElementNodes(loc, *element, bNodes, s);
                        answer.assemble(loc, mat);
                    }
                }
            }
        }
    }

    if ( domain->hasContactManager() ) {
        OOFEM_ERROR("Contant problems temporarily deactivated");
        //domain->giveContactManager()->assembleTangentFromContacts(answer, tStep, type, s, s);
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    answer.assembleBegin();
    answer.assembleEnd();
}


} // end namespace oofem
