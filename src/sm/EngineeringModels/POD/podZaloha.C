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

#include "../sm/EngineeringModels/POD/pod.h"

#include "../sm/EngineeringModels/nlinearstatic.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/structuralelementevaluator.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "error.h"
#include "verbose.h"
#include "sparsenonlinsystemnm.h"
#include "nrsolver.h"
#include "calmls.h"
#include "outputmanager.h"
#include "datastream.h"
#include "classfactory.h"
#include "timer.h"
#include "contextioerr.h"
#include "sparsemtrx.h"
#include "errorestimator.h"
#include "mathfem.h"
#include "dofmanager.h"
#include "dof.h"

#include "util.h"

#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/Micromorphic/micromorphicms.h"
#include "material.h"

#include "gausspoint.h"




#include "engngm.h"
#include "oofemtxtdatareader.h"
#include "skylineu.h"
#include "node.h"
#include "sparsemtrx.h"
#include "exportmodulemanager.h"

namespace oofem {
REGISTER_EngngModel(POD);

  POD :: POD(int i, EngngModel *_master) : NonLinearStatic(i, _master)
{

  equationNumbering = NULL;
  //epsilonPOD = 1.e-8;
 
}


POD :: ~POD()
{
  delete reducedState_dofs;
  delete reducedState_stress;
  if(hyperReductionFlag){
    delete hyperReduction;
  }
  if(plotPodBasisFlag == true){
    delete podVTKXMLExportModule;
  }
  if(equationNumbering) {
    delete equationNumbering;
  }

}


int
POD :: instanciateYourself(DataReader *dr, InputRecord *ir, const char *dataOutputFileName, const char *desc)
{
    int result;
    result = NonLinearStatic :: instanciateYourself(dr, ir, dataOutputFileName, desc);
    ir->finish();
    // instanciate slave problems
    result &= this->instanciateSlaveProblems();
    return result;
}



IRResultType
POD :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    NonLinearStatic :: initializeFrom(ir);

    performAnalysisFlag = true;
    IR_GIVE_OPTIONAL_FIELD(ir, this->performAnalysisFlag, _IFT_POD_performAnalysis);

    performSnapshotsFlag = true;
    IR_GIVE_OPTIONAL_FIELD(ir, this->performSnapshotsFlag, _IFT_POD_performSnapshots);
    if(performSnapshotsFlag == true) {
      IR_GIVE_FIELD(ir, inputStreamNames, _IFT_POD_slaveprob);
    }

    nReducedModes = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->nReducedModes, _IFT_POD_numberOfReducedModes);



    epsOrtho = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->epsOrtho, _IFT_POD_epsortho);
   
    plotPodBasisFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, this->plotPodBasisFlag, _IFT_POD_plotPodBasis);



    if(plotPodBasisFlag == true){
      podVTKXMLExportModule = new PODVTKXMLExportModule(1, this);
      podVTKXMLExportModule->initializeFrom(ir);
    }



    hyperReductionFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, this->hyperReductionFlag, _IFT_POD_hr);

    if(hyperReductionFlag == true){

      hyperReduction = new HyperReduction();
      hyperReduction->initializeFrom(ir);

      this->computeResultsOutsideRIDFlag = false;
      IR_GIVE_OPTIONAL_FIELD(ir, this->computeResultsOutsideRIDFlag, _IFT_POD_computeResultsOutsideRIDFlag);
      if(computeResultsOutsideRIDFlag){
	IR_GIVE_FIELD(ir, this->outOfRIDElementSet, _IFT_POD_outOfRIDElementSet);
      }
	
    }


    evaluateErrorFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, this->evaluateErrorFlag, _IFT_POD_evaluateError);

    


    initReducedStateFromFileFlag = false;
    saveReducedStateToFileFlag = false;


    IR_GIVE_OPTIONAL_FIELD(ir, this->initReducedStateFromFileFlag, _IFT_POD_initReducedStateFromFile);
    IR_GIVE_OPTIONAL_FIELD(ir, this->saveReducedStateToFileFlag, _IFT_POD_saveReducedStateToFile);

    if( performSnapshotsFlag == false && initReducedStateFromFileFlag == false) {
      OOFEM_ERROR("No Snapshots, neither init file");
    }

    separateBasisFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, separateBasisFlag, _IFT_POD_separateBasis);

    numberOfDofGroups = 1;

    if(separateBasisFlag) {
      int index = 0;
      IR_GIVE_FIELD(ir, nDofGroups, _IFT_POD_nDofGroups);
      IR_GIVE_FIELD(ir, dofGroups, _IFT_POD_dofGroups);
      numberOfDofGroups = nDofGroups.giveSize();
      dofGroupArray = new IntArray[numberOfDofGroups];

      for(int i = 0; i < numberOfDofGroups; i ++) {
	dofGroupArray[i].resize(nDofGroups.at(i+1));
	for(int j = 1; j <= nDofGroups.at(i+1); j++) {
	  index++;
	  dofGroupArray[i].at(j) = dofGroups.at(index);
	}
      }     
    }



    reducedState_dofs = new ReducedState[nDofGroups];
    reducedState_stress = new ReducedState[nDofGroups];



    if(initReducedStateFromFileFlag) {


      IR_GIVE_FIELD(ir, dofsReducedStateInputFileName, _IFT_POD_dofsReducedStateInputFileName);
      IR_GIVE_FIELD(ir, stressReducedStateInputFileName, _IFT_POD_stressReducedStateInputFileName);
      /// @todo add separate basis reading of files

      FILE *file = NULL;
      file = fopen(dofsReducedStateInputFileName.c_str(), "rb");
      DataStream *stream = new FileDataStream(file);
      reducedState_dofs->restoreYourself(*stream);
      delete stream;
      file = NULL;
      file = fopen(stressReducedStateInputFileName.c_str(), "rb");
      stream = new FileDataStream(file);
      reducedState_stress->restoreYourself(*stream);
      delete stream;
      
    }

          /// @todo add separate basis writing into files

    if(saveReducedStateToFileFlag) {
 
      IR_GIVE_FIELD(ir, dofsReducedStateOutputFileName, _IFT_POD_dofsReducedStateOutputFileName);
      IR_GIVE_FIELD(ir, stressReducedStateOutputFileName, _IFT_POD_stressReducedStateOutputFileName);
    }


    IR_GIVE_OPTIONAL_FIELD(ir, this->testFlag, _IFT_POD_testFlag);



    
	
    return IRRT_OK;	
    

}




  
void 
POD :: solveYourself()
{
  // compute POD reduced basis
  if(performSnapshotsFlag) {
    this->computeReducedBasis();
  }

  this->buildReducedDomain();
  //print dof basis
  if(plotPodBasisFlag) {
    TimeStep *tStep;
    podVTKXMLExportModule->initialize();
    podVTKXMLExportModule->doOutput(tStep, this->reducedBasisMatrix.giveNumberOfColumns());
  }

  this->printPodBasisMatrix();
 

  StructuralEngngModel :: solveYourself();
}
  
  
void
POD :: solveYourselfAt(TimeStep *tStep)
{
  if(performAnalysisFlag){
    if(testFlag)
      proceedStep2(1,tStep);
    else
      proceedStep(1, tStep);
  }
}
  
  
void
POD :: updateYourself(TimeStep *tStep) 
{
  if(hyperReductionFlag) {
    if(computeResultsOutsideRIDFlag) {
      this->postProcessResults(tStep);
    }
  }
  StructuralEngngModel :: updateYourself(tStep);
}

void
POD :: terminate(TimeStep *tStep)
{

  if(evaluateErrorFlag){
    this->printError(tStep);
  }
  // this should be after postprocessing, but so far leads to an error
  this->printReactionForces(tStep, 1);
  this->doStepOutput(tStep);
  // update load vectors before storing context
  fflush( this->giveOutputStream() );
  this->updateLoadVectors(tStep);
  this->saveStepContext(tStep); 
}

void
POD :: postProcessResults(TimeStep *tStep)
{
  
  // first activate elements in postProcessElementSet
  for (int i = 1; i <= outOfRIDElementSet.giveSize(); i++) {
    int setId = outOfRIDElementSet.at(i);
    IntArray elementNumbers = this->giveDomain(1)->giveSet(setId)->giveElementList();
    for(int j = 1; j <= elementNumbers.giveSize(); j++) {
      this->giveDomain(1)->giveElement(elementNumbers.at(j))->activateYourself();
    }
  }

  // renumber equations
  this->forceEquationNumberingPostProcessing(1);

  IntArray freeDofs(this->reducedState_dofs->giveReducedBasisMatrix().giveNumberOfRows());
  freeDofs.zero();
  int numberOfFreeDofs = 0;
  for(int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++) {
    IntArray dofIDArray;
    this->giveDomain(1)->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
    for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
      if(!this->giveDomain(1)->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	numberOfFreeDofs++;
	freeDofs.at(numberOfFreeDofs) = (i-1)*dofIDArray.giveSize() + k;
      }
    }
  }

  freeDofs.resizeWithValues(numberOfFreeDofs);
  IntArray allColumns(this->reducedState_dofs->giveReducedBasisMatrix().giveNumberOfColumns());

  for(int i = 1; i <= allColumns.giveSize(); i++){
    allColumns.at(i) = i;
  }

  

  rbPostprocess.beSubMatrixOf(this->reducedState_dofs->giveReducedBasisMatrix(),freeDofs, allColumns);
  //  V.beSubMatrixOf(this->reducedState_dofs->giveReducedBasisMatrix(),freeDofs, allColumns);

  totalDisplacement.beProductOf(rbPostprocess,totalReducedCoordinate);
  // compute internal forces, strains and stresses are computed during the step
  FloatArray iF;
  this->giveInternalForces(iF, true, 1, tStep);

}




void
POD :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    if ( ( di == 1 ) && ( tStep == this->giveCurrentStep() ) ) {
      if(hyperReduction && computeResultsOutsideRIDFlag) {
	/*	FloatArray initialLV, incrementalLV;
		initialLV.beProductOf(reducedBasisMatrix,initialLoadVector);
		initialLoadVector.beTProductOf( this->reducedState_dofs->giveReducedBasisMatrix(),initialLV);
		incrementalLV.beProductOf(reducedBasisMatrix,incrementalLoadVector);
		incrementalLoadVector.beProductOf(rbPostprocess,incrementalLV);
		reactions = initialLoadVectorOfPrescribed;
		reactions.add(loadLevel, incrementalLoadVectorOfPrescribed);
	*/

      }	else {
        reactions = initialLoadVectorOfPrescribed;
        reactions.add(loadLevel, incrementalLoadVectorOfPrescribed);
      }
    } else {
        OOFEM_ERROR("unable to respond due to invalid solution step or domain");
    }
}






void
POD :: printError(TimeStep *tStep)
{
  // test if solution step output is active
  if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(tStep) ) {
    return;
  }

  std :: string fileName = giveOutputBaseFileName() + ".err";
  FILE *errorOutFile;
  if ( ( errorOutFile = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str());
    }
  //  fprintf( errorOutFile, "\n\n\t ROM Error \n\n\n");
  //fprintf( errorOutFile, "\tError %.4e\n", rInfNorm );
  fprintf( errorOutFile, "%.4e \n", rInfNorm);
  fprintf( errorOutFile, "%.4e \n", r2Norm);
  fclose( errorOutFile);

}

 
 
void 
POD :: printPodBasisMatrix()  
{

  std :: string fileName = giveOutputBaseFileName() + ".basis";
  FILE *podFileName;
  if ( ( podFileName = fopen(fileName.c_str(), "w") ) == NULL ) {
      OOFEM_ERROR("failed to open file %s", fileName.c_str());
    }
  this->reducedState_dofs->giveReducedBasisMatrix().printYourself(podFileName);
    //  this->reducedBasisMatrix.printYourself(podFileName);
}

  /*
  for(int i =0; i < reducedBasisMatrix.giveSize(); i++){
    std :: string fileName = giveOutputBaseFileName() + ".basis";
    FILE *errorOutFile;
    if ( ( errorOutFile = fopen(fileName.c_str(), "w") ) == NULL ) {
      OOFEM_ERROR("failed to open file %s", fileName.c_str());
    }

  }

  
  //  fprintf( errorOutFile, "\n\n\t ROM Error \n\n\n");
  //fprintf( errorOutFile, "\tError %.4e\n", rInfNorm );
  fprintf( errorOutFile, "%.4e\n", rInfNorm );
  fclose( errorOutFile);


 for ( int j = 1; j <= size; j++ ) {
        DofIDItem id = ( DofIDItem ) dofIDMask.at(j);
	if ( id == Undef ) {
	  answer.at(j) = 0.;
	} else if ( dman->hasDofID(id) ) {
            // primary variable available directly in DOF-manager
	  POD *pod = static_cast<POD*> (emodel);
	  if(pod)
	    answer.at(j) = pod->giveReducedBasisValue(tStep,dman, dofIDMask.at(j), iSnapshot);
//reducedBasisMatrix.at(,1);
	}
    }

}
 
  */





void
POD :: proceedStep(int di, TimeStep *tStep)
{



 if ( initFlag ) {
        //
        // first step  create space for stiffness Matrix
        //
        if ( !stiffnessMatrix ) {
            stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        }

        if ( !stiffnessMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        if ( nonlocalStiffnessFlag ) {
            if ( !stiffnessMatrix->isAsymmetric() ) {
                OOFEM_ERROR("stiffnessMatrix does not support asymmetric storage");
            }
        }

        stiffnessMatrix->buildInternalStructure( this, di, EModelDefaultEquationNumbering() );
    }

#if 0
    if ( ( mstep->giveFirstStepNumber() == tStep->giveNumber() ) ) {
 #ifdef VERBOSE
        OOFEM_LOG_INFO("Resetting load level\n");
 #endif
        if ( mstepCumulateLoadLevelFlag ) {
            cumulatedLoadLevel += loadLevel;
        } else {
            cumulatedLoadLevel = 0.0;
        }
        this->loadLevel = 0.0;
    }
#endif

    if ( loadInitFlag || controlMode == nls_directControl ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling reference load\n");
#endif
        //
        // assemble the incremental reference load vector
        //

	FloatArray iLV, iLVP;

        this->assembleIncrementalReferenceLoadVectors(iLV,incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, this->giveDomain(di), tStep);

	incrementalLoadVector.beTProductOf(reducedBasisMatrix, iLV);

        loadInitFlag = 0;
    }

    if ( tStep->giveNumber() == 1 ) {
        int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();

	totalReducedCoordinate.resize(reducedBasisMatrix.giveNumberOfColumns());
	totalReducedCoordinate.zero();
	incrementOfReducedCoordinate.resize(reducedBasisMatrix.giveNumberOfColumns());
	incrementOfReducedCoordinate.zero();
    }

    //
    //    ->   BEGINNING OF LOAD (OR DISPLACEMENT) STEP  <-
    //

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    //
    // call numerical model to solve arise problem
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving       [step number %5d.%d, time = %e]\n\n", tStep->giveNumber(), tStep->giveVersion(), tStep->giveIntrinsicTime() );
#endif

    FloatArray extrapolatedForces;
    FloatArray *extrapolatedForcesPtr = &extrapolatedForces;
    if ( this->initialGuessType == IG_Original ) {
      incrementOfDisplacement.zero();
      this->assembleExtrapolatedForces( extrapolatedForces, tStep, ElasticStiffnessMatrix, this->giveDomain(di) );
      extrapolatedForces.negated();
    } else {
      incrementOfDisplacement.zero();
      extrapolatedForcesPtr = NULL;
    }


    
    



    //totalDisplacement.printYourself();
    if ( initialLoadVector.isNotEmpty() ) {
      numMetStatus = nMethod->solveHR(*stiffnessMatrix, incrementalLoadVector, & initialLoadVector, extrapolatedForcesPtr,
                                      totalReducedCoordinate, incrementOfReducedCoordinate, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep, reducedBasisMatrix);
    } else {
      numMetStatus = nMethod->solveHR(*stiffnessMatrix, incrementalLoadVector, NULL, extrapolatedForcesPtr,
                                      totalReducedCoordinate, incrementOfReducedCoordinate, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep, reducedBasisMatrix);
    }
    ///@todo Martin: ta bort!!!
    //this->updateComponent(tStep, NonLinearLhs, this->giveDomain(di));
    ///@todo Use temporary variables. updateYourself() should set the final values, while proceedStep should be callable multiple times for each step (if necessary). / Mikael
    OOFEM_LOG_RELEVANT("Equilibrium reached at load level = %f in %d iterations\n", cumulatedLoadLevel + loadLevel, currentIterations);
    prevStepLength =  currentStepLength;
    

}




void
POD :: proceedStep2(int di, TimeStep *tStep)
{

 if ( initFlag ) {
        //
        // first step  create space for stiffness Matrix
        //
        if ( !stiffnessMatrix ) {
            stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        }

        if ( !stiffnessMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        if ( nonlocalStiffnessFlag ) {
            if ( !stiffnessMatrix->isAsymmetric() ) {
                OOFEM_ERROR("stiffnessMatrix does not support asymmetric storage");
            }
        }

        stiffnessMatrix->buildInternalStructure( this, di, *equationNumbering);
    }

#if 0
    if ( ( mstep->giveFirstStepNumber() == tStep->giveNumber() ) ) {
 #ifdef VERBOSE
        OOFEM_LOG_INFO("Resetting load level\n");
 #endif
        if ( mstepCumulateLoadLevelFlag ) {
            cumulatedLoadLevel += loadLevel;
        } else {
            cumulatedLoadLevel = 0.0;
        }
        this->loadLevel = 0.0;
    }
#endif

    if ( loadInitFlag || controlMode == nls_directControl ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling reference load\n");
#endif
        //
        // assemble the incremental reference load vector
        //

	FloatArray iLV, iLVP;

        this->assembleIncrementalReferenceLoadVectors2(iLV,incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, this->giveDomain(di), tStep);

	incrementalLoadVector.beTProductOf(reducedBasisMatrix, iLV);

        loadInitFlag = 0;
    }

    if ( tStep->giveNumber() == 1 ) {
        int neq = this->giveNumberOfDomainEquations( 1, *equationNumbering);
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();

	totalReducedCoordinate.resize(reducedBasisMatrix.giveNumberOfColumns());
	totalReducedCoordinate.zero();
	incrementOfReducedCoordinate.resize(reducedBasisMatrix.giveNumberOfColumns());
	incrementOfReducedCoordinate.zero();
    }

    //
    //    ->   BEGINNING OF LOAD (OR DISPLACEMENT) STEP  <-
    //

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    //
    // call numerical model to solve arise problem
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving       [step number %5d.%d, time = %e]\n\n", tStep->giveNumber(), tStep->giveVersion(), tStep->giveIntrinsicTime() );
#endif

    FloatArray extrapolatedForces;
    FloatArray *extrapolatedForcesPtr = &extrapolatedForces;
    if ( this->initialGuessType == IG_Original ) {
      incrementOfDisplacement.zero();
      this->assembleExtrapolatedForces( extrapolatedForces, tStep, ElasticStiffnessMatrix, this->giveDomain(di) );
      extrapolatedForces.negated();
    } else {
      incrementOfDisplacement.zero();
      extrapolatedForcesPtr = NULL;
    }


    
    



    //totalDisplacement.printYourself();
    if ( initialLoadVector.isNotEmpty() ) {
      numMetStatus = nMethod->solveHR(*stiffnessMatrix, incrementalLoadVector, & initialLoadVector, extrapolatedForcesPtr,
                                      totalReducedCoordinate, incrementOfReducedCoordinate, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep, reducedBasisMatrix);
    } else {
      numMetStatus = nMethod->solveHR(*stiffnessMatrix, incrementalLoadVector, NULL, extrapolatedForcesPtr,
                                      totalReducedCoordinate, incrementOfReducedCoordinate, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep, reducedBasisMatrix);
    }
    ///@todo Martin: ta bort!!!
    //this->updateComponent(tStep, NonLinearLhs, this->giveDomain(di));
    ///@todo Use temporary variables. updateYourself() should set the final values, while proceedStep should be callable multiple times for each step (if necessary). / Mikael
    OOFEM_LOG_RELEVANT("Equilibrium reached at load level = %f in %d iterations\n", cumulatedLoadLevel + loadLevel, currentIterations);
    prevStepLength =  currentStepLength;
    

}



 

void
POD :: assembleIncrementalReferenceLoadVectors2(FloatArray &_incrementalLoadVector,
                                                           FloatArray &_incrementalLoadVectorOfPrescribed,
                                                           SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                           Domain *sourceDomain, TimeStep *tStep)
{
  _incrementalLoadVector.resize( sourceDomain->giveEngngModel()->giveNumberOfDomainEquations( sourceDomain->giveNumber(), *equationNumbering ));
    _incrementalLoadVector.zero();
    _incrementalLoadVectorOfPrescribed.resize( sourceDomain->giveEngngModel()->giveNumberOfDomainEquations( sourceDomain->giveNumber(), EModelDefaultPrescribedEquationNumbering() ) );
    _incrementalLoadVectorOfPrescribed.zero();

    if ( _refMode == SparseNonLinearSystemNM :: rlm_incremental ) {
        ///@todo This was almost definitely wrong before. It never seems to be used. Is this code even relevant?
        this->assembleVector(_incrementalLoadVector, tStep, ExternalForceAssembler(),
                             VM_Incremental, *equationNumbering, sourceDomain);

        this->assembleVector(_incrementalLoadVectorOfPrescribed, tStep, ExternalForceAssembler(),
                             VM_Incremental, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    } else {
        this->assembleVector(_incrementalLoadVector, tStep, ExternalForceAssembler(),
                             VM_Total, *equationNumbering, sourceDomain);

        this->assembleVector(_incrementalLoadVectorOfPrescribed, tStep, ExternalForceAssembler(),
                             VM_Total, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    }

    this->updateSharedDofManagers(_incrementalLoadVector, *equationNumbering, LoadExchangeTag);
}


int
POD :: forceEquationNumbering2(int id)
{
    // forces equation renumbering for current time step
    // intended mainly for problems with changes of static system
    // during solution
    // OUTPUT:
    // sets this->numberOfEquations and this->numberOfPrescribedEquations and returns this value

    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    for ( auto &node : domain->giveDofManagers() ) {
      if(hyperReduction->giveSelectedNodes().at(node->giveNumber()) == 1) {
	node->askNewEquationNumbers(currStep);
      }
    }
    
    for ( auto &elem : domain->giveElements() ) {
      if(elem->isActivated()) {
	int nnodes = elem->giveNumberOfInternalDofManagers();
	for ( int k = 1; k <= nnodes; k++ ) {
	  elem->giveInternalDofManager(k)->askNewEquationNumbers(currStep);
	}
      }
    }
    
     equationNumberingCompleted = 1;
    
    return domainNeqs.at(id);
}

int
POD :: forceEquationNumberingPostProcessing(int id)
{
    // forces equation renumbering for current time step
    // intended mainly for problems with changes of static system
    // during solution
    // OUTPUT:
    // sets this->numberOfEquations and this->numberOfPrescribedEquations and returns this value

    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    for ( auto &node : domain->giveDofManagers() ) {
      node->askNewEquationNumbers(currStep);
    }
    
    for ( auto &elem : domain->giveElements() ) {
   	int nnodes = elem->giveNumberOfInternalDofManagers();
	for ( int k = 1; k <= nnodes; k++ ) {
	  elem->giveInternalDofManager(k)->askNewEquationNumbers(currStep);
	}
    }
    
    
     equationNumberingCompleted = 1;
    
    return domainNeqs.at(id);
}

void
POD :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
//
// updates some component, which is used by numerical method
// to newly reached state. used mainly by numerical method
// when new tangent stiffness is needed during finding
// of new equilibrium stage.
//
{
    switch ( cmpn ) {
    case NonLinearLhs:
      {
	if ( stiffMode == nls_tangentStiffness ) {
	  stiffnessMatrix->zero(); // zero stiffness matrix
#ifdef VERBOSE
	  OOFEM_LOG_DEBUG("Assembling tangent stiffness matrix\n");
#endif
	  stiffnessMatrix->buildInternalStructure( this, d->giveNumber(), EModelDefaultEquationNumbering() );
	  this->assemble(*stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
			 EModelDefaultEquationNumbering(), d);
	  
	} else if ( ( stiffMode == nls_secantStiffness ) || ( stiffMode == nls_secantInitialStiffness && initFlag ) ) {
#ifdef VERBOSE
	  OOFEM_LOG_DEBUG("Assembling secant stiffness matrix\n");
#endif
	  stiffnessMatrix->buildInternalStructure( this, d->giveNumber(), EModelDefaultEquationNumbering() );
	  stiffnessMatrix->zero(); // zero stiffness matrix
	  this->assemble(*stiffnessMatrix, tStep, TangentAssembler(SecantStiffness),
			 EModelDefaultEquationNumbering(), d);
	  initFlag = 0;
	} else if ( ( stiffMode == nls_elasticStiffness ) && ( initFlag ||
							       ( this->giveMetaStep( tStep->giveMetaStepNumber() )->giveFirstStepNumber() == tStep->giveNumber() ) || (updateElasticStiffnessFlag) ) ) {
#ifdef VERBOSE
	  OOFEM_LOG_DEBUG("Assembling elastic stiffness matrix\n");
#endif
	  stiffnessMatrix->zero(); // zero stiffness matrix
	  stiffnessMatrix->buildInternalStructure( this, d->giveNumber(), EModelDefaultEquationNumbering() );
	  this->assemble(*stiffnessMatrix, tStep, TangentAssembler(ElasticStiffness),
			 EModelDefaultEquationNumbering(), d);
	  initFlag = 0;
	} else {
	// currently no action , this method is mainly intended to
	// assemble new tangent stiffness after each iteration
	// when secantStiffMode is on, we use the same stiffness
	// during iteration process
	}


	
	FloatMatrix KA, ATKA;
	//	SkylineUnsym *skyReducedStiffnessMatrix;
	std :: unique_ptr<SkylineUnsym> srStiffnessMatrix;
	std :: unique_ptr<SkylineUnsym> testMtrx;
	std :: unique_ptr< SparseMtrx > reducedStiffnessMatrix;





	srStiffnessMatrix =  std::unique_ptr<SkylineUnsym>(static_cast<SkylineUnsym*> (stiffnessMatrix.release()));

	
	srStiffnessMatrix->times(reducedBasisMatrix, KA);


	if(hyperReductionFlag) {
	  IntArray interfaceDofs;
	  interfaceDofs = hyperReduction->giveInterfaceDofs();
	  for(int j = 1; j <= reducedBasisMatrix.giveNumberOfColumns(); j++) {	  
	    for(int i = 1; i <= interfaceDofs.giveSize(); i++) {	   
	      KA.at(interfaceDofs.at(i),j) = 0;
	    }
	  }
	}

	ATKA.beTProductOf(reducedBasisMatrix,KA);
	//transform to skyline and to unique_pntr
	srStiffnessMatrix->initializeFromFloatMatrix(ATKA);
	stiffnessMatrix =  std::unique_ptr<SparseMtrx>(static_cast<SparseMtrx*> (srStiffnessMatrix.release()));


      }
      break;
    case InternalRhs:
      {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Updating internal forces\n");
#endif
	FloatArray iF;
	totalDisplacement.beProductOf(reducedBasisMatrix,totalReducedCoordinate );
	incrementOfDisplacement.beProductOf(reducedBasisMatrix,incrementOfReducedCoordinate );
        // update internalForces and internalForcesEBENorm concurrently
        this->giveInternalForces(iF, true, d->giveNumber(), tStep);
	if(hyperReductionFlag) {
	  IntArray interfaceDofs = hyperReduction->giveInterfaceDofs();
	  for(int i = 1; i <= interfaceDofs.giveSize(); i++) {
	    iF.at(interfaceDofs.at(i)) = 0;
	  }
	}
	
	internalForces.beTProductOf(reducedBasisMatrix,iF);
	if(evaluateErrorFlag) {
	  rInfNorm = 0;
	  r2Norm = 0;
	  for (int i = 1; i<= iF.giveSize(); i++) {
	    r2Norm += iF.at(i) * iF.at(i);
	    if(fabs(iF.at(i)) > rInfNorm ) {
	      rInfNorm = fabs(iF.at(i));
	     
	    }
	  }
	  r2Norm = sqrt(r2Norm);
	}



        break;
      }
    default:
        OOFEM_ERROR("Unknown Type of component.");
    }

}




void
POD :: computeReducedBasis()
{

  // Solve slave problems and save their solution as snapshots 
  int nSnapshots = 0;
  int smstep = 1, sjstep = 1;
  EngngModel *sp;
  nStressComponents = 0;
  for ( int i = 1; i <= this->giveNumberOfSlaveProblems(); i++ ) {
    sp = this -> giveSlaveProblem(i);
    FILE *out = sp->giveOutputStream();

    //	sp -> solveYourself();
    sp->giveTimer()->startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    if ( sp->giveCurrentStep() ) {
      smstep = sp->giveCurrentStep()->giveMetaStepNumber();
      sjstep = sp->giveMetaStep(smstep)->giveStepRelativeNumber( sp->giveCurrentStep()->giveNumber() ) + 1;
    } 
    // solve slave problems to get snapshots 
    for ( int imstep = smstep; imstep <= nMetaSteps; imstep++, sjstep = 1 ) { //loop over meta steps
      MetaStep *activeMStep = sp->giveMetaStep(imstep);
      // update state according to new meta step
      sp->initMetaStepAttributes(activeMStep);
      int nTimeSteps = activeMStep->giveNumberOfSteps();
      for ( int jstep = sjstep; jstep <= nTimeSteps; jstep++ ) { //loop over time steps
	nSnapshots++;
	sp->giveTimer()->startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
	sp->giveTimer()->initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
	sp->giveNextStep();
	// renumber equations if necessary. Ensure to call forceEquationNumbering() for staggered problems
	if ( sp->requiresEquationRenumbering( sp->giveCurrentStep() ) ) {
	  sp->forceEquationNumbering();
	}
	OOFEM_LOG_DEBUG("Number of equations %d\n", sp->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering()) );

	sp->solveYourselfAt( sp->giveCurrentStep() );
	sp->updateYourself( sp->giveCurrentStep() );

	sp->giveTimer()->stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
	sp->terminate( sp->giveCurrentStep() );
	
	double _steptime = sp->giveSolutionStepTime();
	OOFEM_LOG_INFO("EngngModel info: user time consumed by solution step %d: %.2fs\n", sp->giveCurrentStep()->giveNumber(), _steptime);
	fprintf(out, "\nUser time consumed by solution step %d: %.3f [s]\n\n",sp->giveCurrentStep()->giveNumber(), _steptime);
	

	FloatArray FEM_vars;
	// take the snapshot for dofs	
	this->takeSnapshot_dofs(FEM_vars, sp->giveCurrentStep(), sp->giveDomain(1), dofIDMatrix );
	// add snapshot to dof basis
	//	reducedState_dofs->subspaceExpansion_dofs(FEM_vars);
	reducedState_dofs->subspaceExpansion(FEM_vars);
	// take the snapshot for stresses
	this->takeSnapshot_stress(FEM_vars, sp->giveCurrentStep(), sp->giveDomain(1), nStressComponents);
	// add snapshot to dof basis
	//	reducedState_stress->subspaceExpansion_stress(FEM_vars);
	reducedState_stress->subspaceExpansion(FEM_vars);
	reducedState_stress->setNumberOfStressComponents(nStressComponents);

      }
    }
  }
  // subspace selection for dofs and stresses
  //  bool ret = reducedState_dofs->subspaceSelection_dofs();
  bool ret = reducedState_dofs->subspaceSelection();
  if(ret)
    //    ret = reducedState_stress->subspaceSelection_stress();
    ret = reducedState_stress->subspaceSelection();
  if(!ret)
    OOFEM_ERROR("Huhu error in subspace selection");

 

  
  if(saveReducedStateToFileFlag) {
    // store reduced basis into file
    FILE *file = NULL;
    file = fopen(dofsReducedStateOutputFileName.c_str(), "wb");
    //    file = fopen(dofsReducedStateOutputFileName.c_str(), "rs.out.dofs");
    DataStream *stream = new FileDataStream(file);
    reducedState_dofs->storeYourself(*stream);
    delete stream;
    
    file = NULL;
    file = fopen(stressReducedStateOutputFileName.c_str(), "wb");
    //    file = fopen(dofsReducedStateOutputFileName.c_str(), "rs.out.stress");
    stream = new FileDataStream(file);
    reducedState_stress->storeYourself(*stream);
    delete stream;
  }

}


void 
POD :: buildReducedDomain()
{
  // build reduced integration domain
  if(hyperReductionFlag) {    
    hyperReduction->initializeYourself(this->giveDomain(1));
    // find reduced domain nodes
    hyperReduction->computeReducedDomainNodes_dofs(this->giveCurrentStep(), reducedState_dofs->giveReducedBasisMatrix());
    hyperReduction->computeReducedDomainNodes_stress(this->giveCurrentStep(), reducedState_stress->giveReducedBasisMatrix(),reducedState_stress->giveNumberOfStressComponents());
    // build reduced domain
    //@todo uncomment following line
    if(testFlag)
      hyperReduction->buildReducedDomain2(this->giveCurrentStep(),reducedState_dofs->giveReducedBasisMatrix());
    else
      hyperReduction->buildReducedDomain(this->giveCurrentStep(),reducedState_dofs->giveReducedBasisMatrix());
    // get hyperReduced basis matrix
    reducedBasisMatrix = hyperReduction->giveHyperReducedBasisMatrix();

    if(nReducedModes != 0) {
      int nDofs = reducedBasisMatrix.giveNumberOfRows();
      int nModes = reducedBasisMatrix.giveNumberOfColumns();
      if(nModes >= nReducedModes) {
	reducedBasisMatrix.resize(nDofs,nReducedModes);      
      }
    }
    
    // set reduced domain as the engngm domain
    //@todo uncomment following lines
    if(!testFlag) {
      this->setDomain(1,hyperReduction->giveReducedDomain(),0);
      // reinitialize export module    
      exportModuleManager->reInitialize();
      // renumber force equations - is it necessary?
      this->forceEquationNumbering(1);
    }
    if(testFlag) {
      equationNumbering = new ReducedDomainNumberingScheme(hyperReduction->giveSelectedNodes());
      equationNumbering->init(this->giveDomain(1), this->giveCurrentStep());
      this->forceEquationNumbering2(1);
    }
  } else {
    reducedBasisMatrix = reducedState_dofs->giveReducedBasisMatrix();
  }

}


void
POD :: takeSnapshot_dofs(FloatArray &answer, TimeStep *tStep, Domain *domain, std:: vector<IntArray> &dofIDMatrix )
{

  IntArray dofIDArray;
  FloatArray vec;
  answer.resize(0);

  int nDofs = 0;
  for (int iNode = 1; iNode <= domain->giveNumberOfDofManagers(); iNode ++ ) {	
    DofManager* dofManager = domain->giveDofManager(iNode);
    dofManager->giveCompleteMasterDofIDArray(dofIDArray);
    
    nDofs += dofIDArray.giveSize();
    IntArray ia;
    ia.followedBy(iNode);
    ia.followedBy(dofIDArray);
    ia.followedBy(nDofs);
    
    dofIDMatrix.push_back(ia);
    dofManager->giveUnknownVector(vec, dofIDArray, VM_Total, tStep, true);
    answer.append(vec);
  }
}

void
POD :: takeSnapshot_stress(FloatArray &answer, TimeStep *tStep, Domain *domain, int &stressSize)
{

  IntArray dofIDArray;
  FloatArray vec;
  answer.resize(0);

  for (int iElement = 1; iElement <= domain->giveNumberOfElements(); iElement++ ) {
    Element *element = (domain->giveElement(iElement));
    /*StructuralElementEvaluator *see;
    if(!se){
      see = dynamic_cast< StructuralElementEvaluator* >(domain->giveElement(iElement));
      if(!see)
	OOFEM_ERROR("POD works only for structural problems")
    } 
    
    */



    for( GaussPoint *gp: *element->giveDefaultIntegrationRulePtr()) {
      
      Material *m = gp->giveMaterial();
      MaterialStatus *ms = m->giveStatus(gp);
      //->giveMaterialStatus();
      
      MicromorphicMaterialStatus *mms = dynamic_cast< MicromorphicMaterialStatus * >(ms); 
      StructuralMaterialStatus *sms = dynamic_cast< StructuralMaterialStatus * >(ms); 

      FloatArray stress(0);
      if(mms) {
	FloatArray sigma, M;
	sigma = mms -> giveStressVector();
	M = mms -> giveMicromorphicStressGrad();
	stress.append(sigma);
	stress.append(M);
      } else if(sms) {
	stress = sms->giveStressVector();
      }
      

      stressSize = stress.giveSize();
      answer.append(stress);

    }
    
  }
}



EngngModel *
POD :: giveSlaveProblem(int i)
{
    if ( ( i > 0 ) && ( i <= this->giveNumberOfSlaveProblems() ) ) {
        return this->emodelList[i-1].get();
    } else {
        OOFEM_ERROR("Undefined problem");
    }

    return NULL;
}


int
POD :: instanciateSlaveProblems()
{
    //first instantiate master problem if defined
    emodelList.resize(inputStreamNames.size());
  
    for ( int i = 1; i <= (int)inputStreamNames.size(); i++ ) {
       
        OOFEMTXTDataReader dr( inputStreamNames [ i - 1 ] );
	std :: unique_ptr< EngngModel >prob( InstanciateProblem(& dr, this->pMode, this->contextOutputMode, NULL) );
        emodelList[i-1] = std::move(prob);
    }

    return 1;
}

 

double 
POD :: giveReducedBasisValue(TimeStep *tStep, DofManager *dofManager, int dofID, int iSnapshot)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{

  // numbering from 0
  int nodeNumber = dofManager->giveNumber()-1;
  int iDof = 0;
  for (int i = 2; i < (dofIDMatrix.at(nodeNumber)).giveSize();i++) {
    if((dofIDMatrix.at(nodeNumber)).at(i) == dofID) {
      iDof = i-1;
      break;
    }
  }
  
  if(iDof == 0)
    return 0;
  else {
    int shift = 0;
    if(nodeNumber != 0) {
      shift = (dofIDMatrix.at(nodeNumber-1)).at((dofIDMatrix.at(nodeNumber-1)).giveSize());
    }

    int index = iDof + shift;
    //    if(iSnapshot > reducedState_dofs->giveReducedBasisMatrix_dofs().giveNumberOfColumns()){
    if(iSnapshot > reducedState_dofs->giveReducedBasisMatrix().giveNumberOfColumns()){
      return 0;
    } else{      
      //      return reducedState_dofs->giveReducedBasisMatrix_dofs().at(index, iSnapshot);
      return reducedState_dofs->giveReducedBasisMatrix().at(index, iSnapshot);
    }
  }
 
}


double 
POD :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
 
  int eq = dof->__giveEquationNumber();
  if ( eq == 0 ) {
    return 0.;
  }
  


    if ( tStep != this->giveCurrentStep() ) {
        OOFEM_ERROR("unknown time step encountered");
        return 0.;
    }

    switch ( mode ) {
    case VM_Incremental:
        // return incrementOfDisplacement -> at(eq);
        // return nMethod-> giveUnknownComponent(IncrementOfSolution, eq);
        if ( incrementOfDisplacement.isNotEmpty() ) {
            return incrementOfDisplacement.at(eq);
        } else {
            return 0.;
        }

    case VM_Total:
        if ( totalDisplacement.isNotEmpty() ) {
            return totalDisplacement.at(eq);
        } else {
            return 0.;
        }

    case VM_Velocity:
        if ( incrementOfDisplacement.isNotEmpty() ) {
            return incrementOfDisplacement.at(eq) / tStep->giveTimeIncrement();
        } else {
            return 0.;
        }

    default:
        OOFEM_ERROR("Unknown is of undefined ValueModeType for this problem");
    }

    return 0.0;
}





void 
POD :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
  if(hyperReductionFlag && !computeResultsOutsideRIDFlag) {
    DofManager *dMan = iDof -> giveDofManager();
    if(hyperReduction->giveSelectedNodes().at(dMan->giveNumber())){
         iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total);
    }
  } else {
   iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total);
 }
}



} // end namespace oofem
