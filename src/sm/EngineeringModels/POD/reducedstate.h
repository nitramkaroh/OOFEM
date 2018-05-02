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

#ifndef reducedstate_h
#define reducedstate_h

#include "inputrecord.h"
#include "floatmatrix.h"

#include <string>


#define _IFT_ReducedState_Name "reducedstate"
#define _IFT_ReducedState_tolortho "tolortho"



namespace oofem {
/**
 * huhu popis ReducedState class
**/
class ReducedState 
{
protected:
  /// number of stress components in the problem
  int numberOfStressComponents;
  double epsilonPOD;
  /// reduced basis coordinates
  FloatMatrix rbCoords;
  /*  FloatMatrix rbCoords_dofs;
      FloatMatrix rbCoords_stress;*/
  /// matrix containing reducedBasis
  FloatMatrix reducedBasisMatrix;
  /*  FloatMatrix reducedBasisMatrix_dofs;
  FloatMatrix reducedBasisMatrix_stress;
  */
  /// covariance matrix
  FloatMatrix covarianceMatrix;
  /*  FloatMatrix covarianceMatrix_dofs;
  FloatMatrix covarianceMatrix_stress;
  */
  

  /// input file of reduced state for dofs 
  std :: string dofsReducedStateInputFileName;
  /// input file of reduced state for stress
  std :: string stressReducedStateInputFileName;
  /// output file of reduced state for dofs
  std :: string dofsReducedStateOutputFileName;
  /// output file of reduced state for stress
  std :: string stressReducedStateOutputFileName;


public:
  ReducedState(){;}
  virtual ~ReducedState(){;}

    virtual IRResultType initializeFrom(InputRecord *ir);  
    // identification
    virtual const char *giveClassName() const { return "ReducedState"; }
    virtual const char *giveInputRecordName() const { return _IFT_ReducedState_Name; }
    /*    const FloatMatrix &giveReducedBasisMatrix_dofs() const {return reducedBasisMatrix_dofs;}
  const FloatMatrix &giveReducedBasisMatrix_stress() const {return reducedBasisMatrix_stress;}
    */

    const FloatMatrix &giveReducedBasisMatrix() const {return reducedBasisMatrix;}

  void subspaceExpansion(FloatArray &FEM_vars, double weight = 1) ;
  void subspaceExpansion_dofWeights(FloatArray &FEM_vars, FloatArray &dofWeights, double weight = 1); 
  bool subspaceSelection();
  void storeYourself(DataStream &stream);
  void restoreYourself(DataStream &stream);


  void setNumberOfStressComponents(int nSCp){numberOfStressComponents = nSCp;}

  int giveNumberOfStressComponents() {return numberOfStressComponents;}

  // void subspaceExpansion_dofs(FloatArray &FEM_vars, double weight = 1.0);
  // void subspaceExpansion_stress(FloatArray &FEM_vars, double weight = 1.0);
  // bool subspaceSelection_dofs();
  // bool subspaceSelection_stress();

  /*  void setDofsReducedBasisInputFile(std :: string file){dofsReducedStateInputFileName = file;}
    void setStressReducedBasisInputFile(std :: string file){stressReducedStateInputFileName = file;}
    void setDofsReducedBasisOutputFile(std :: string file){dofsReducedStateOutputFileName = file;}
    void setStressReducedBasisOutputFile(std :: string file){stressReducedStateOutputFileName = file;}
  */

 

protected:
    // void subspaceExpansion(FloatArray &FEM_vars, FloatMatrix &rbM, FloatMatrix &rbCoords, FloatMatrix &cM, double weight) ;
    // bool subspaceSelection(FloatMatrix &rbM, FloatMatrix &rbCoords, FloatMatrix &cM);
    void orthogonalize(FloatArray &answer, const FloatMatrix& A);
    void orthogonalize(FloatArray &answer, const FloatMatrix& A, const FloatArray &dofWeights);

};
} // end namespace oofem
#endif // reducedstate_h
