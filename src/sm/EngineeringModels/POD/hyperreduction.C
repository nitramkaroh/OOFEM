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

#include "../sm/EngineeringModels/POD/hyperreduction.h"
#include "dof.h"
#include "node.h"
#include "dofmanager.h"
#include "element.h"
#include "mathfem.h"
#include "set.h"

#include "iga/iga.h"


namespace oofem {

  HyperReduction :: HyperReduction(Domain *domain): ridNodeSetNumbers(), ridElement2ReduceSetNumbers()
{
  this->domain = domain;
  selectedNodes(0);
  selectedNodes.zero();
  selectedElements(0);
  selectedElements.zero();
}

  HyperReduction :: ~HyperReduction()
{
    delete dof2NodeMap;
}


IRResultType
HyperReduction :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
  
    dof2NodeMap = NULL;

    ridExtension = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, this->ridExtension, _IFT_HyperReduction_ridExtension);        
    IR_GIVE_OPTIONAL_FIELD(ir, this->ridNodeSetNumbers, _IFT_HyperReduction_ridNodeSet);        
    IR_GIVE_OPTIONAL_FIELD(ir, this->ridElement2ReduceSetNumbers, _IFT_HyperReduction_ridElement2ReduceSet);  

    igaModeFlag = false;

    IR_GIVE_OPTIONAL_FIELD(ir, this->igaModeFlag, _IFT_HyperReduction_igaMode);       


    eimDofsSelectionFlag = true;
    IR_GIVE_OPTIONAL_FIELD(ir, this->eimDofsSelectionFlag, _IFT_HyperReduction_EIMDofsSelection); 
    eimStressSelectionFlag = true;
    IR_GIVE_OPTIONAL_FIELD(ir, this->eimStressSelectionFlag, _IFT_HyperReduction_EIMStressSelection);       

    return IRRT_OK;
}
  

void
HyperReduction :: initializeYourself(Domain *d)
{
  this->domain = d;  // initialize selected nodes and elements 
  selectedNodes.resize(d->giveNumberOfDofManagers());
  selectedNodes.zero(); 
  selectedElements.resize(d->giveNumberOfElements());
  selectedElements.zero();



}
   



void
HyperReduction :: computeReducedDomainNodes_dofs(TimeStep* tStep, const FloatMatrix &reducedBasisMatrix) 
{

  int nNodes = domain->giveNumberOfDofManagers();

  IntArray nodesMultiplicity(nNodes);
  nodesMultiplicity.zero();
  this->dof2NodeMap = new std :: map< int, int >();

  int nDofs = 0;
  /// create dof to node map
  for ( int iNode = 1; iNode <= nNodes; iNode++ ) {
        DofManager *dman = domain->giveDofManager(iNode);
        for ( Dof *dof : *dman ) {//= d->giveDofManager(i)->giveDof(j);
	  nDofs++;
	  ( * this->dof2NodeMap ) [ nDofs ] = iNode;
	}
  }    
   



  
  IntArray tmp;
  // perform EIM for dof   
  if(eimDofsSelectionFlag) {
    bool ret = this->EIM(tmp, reducedBasisMatrix);
    if (ret == false) {
      OOFEM_ERROR("solver failed in dof EIM");
    }
  }


  // add Dof_EIM selected nodes to RID
  for (int i = 1; i <= tmp.giveSize(); i++) {
    //@todo this is working if the u is stored first for all nodes then v, and so on
    //    int index = tmp.at(i)%domain->giveNumberOfDofManagers();
    // @todo this is working if all nodes have the same number of dofs
    // @todo should there be +1 or not???
    // @todo uncomment this ???
    //    int index = tmp.at(i)/nDofs;//+1;
    int index = dof2NodeMap->find(tmp.at(i)-1)->second;
    //    int index = dof2NodeMap.at(tmp.at(i));
    //@todo uncomment this
    if(index > selectedNodes.giveSize() || index < 0)
      int ahoj = 1;
    selectedNodes.at(index) = 1;
  }


  /* this nodes are added using ridNodeSet 
  // add nodes with nonzero dirichlet boundary condition
  for(int iNode = 1; iNode <= domain->giveNumberOfDofManagers();iNode++) {
    DofManager *dofManager = domain->giveDofManager(iNode);
    FloatArray vec;
    IntArray dofIDArray;
    dofManager->giveCompleteMasterDofIDArray(dofIDArray);
    dofManager->givePrescribedUnknownVector(vec, dofIDArray,VM_Total, tStep);
    if(vec.giveSize()) {
      for(int i = 1; i<= vec.giveSize(); i++) {
	if(vec.at(i) != 0) {
	  selectedNodes.at(iNode) = 1;
	  break;
	}
      }
    }    
  }
  */

  // this was moved behind rid extension, so that there are not lot of elements around boundary

    for(int i = 1; i <= ridNodeSetNumbers.giveSize(); i++) {
      int setId = ridNodeSetNumbers.at(i);
      IntArray nodesNumbers = this->domain->giveSet(setId)->giveNodeList();
      for(int j = 1; j <= nodesNumbers.giveSize(); j++) {
	selectedNodes.at(nodesNumbers.at(j)) = 1;
      }
    }

}

void
HyperReduction :: computeReducedDomainNodes_stress(TimeStep* tStep, const FloatMatrix &reducedBasisMatrix, int nStressComponents) 
{

  // build indexes gp_to_elem and elem_to_gp
   int count = 0;
   int nElements = this->domain->giveNumberOfElements();
   for (int i = 1; i <= nElements; i++) { 
     count += this->domain->giveElement(i)->giveNumberOfGaussPoints();
       //mesh.elements[i]->num_gp;
   }

   IntArray gp2elem, elem2gp;
   gp2elem.resize(count);
   elem2gp.resize(nElements);

   //ARRAY<int> gp_to_elem(count);
   //elem_to_gp.resize(mesh.elements.size()+1);
   int countGP = 0;
   for (int i = 1; i <= nElements; i++) {
     elem2gp.at(i) = countGP;
     int nGP = this->domain->giveElement(i)->giveNumberOfGaussPoints();
     for (int k = 1; k<= nGP; k++) { 
       gp2elem.at(countGP + k) = i-1;
     }
     countGP += nGP;
   }
   // elem2gp.at(nElements) = countGP;
  
  IntArray tmp;
  // perform EIM for dof   
  if(eimStressSelectionFlag) {
    bool ret = this->EIM(tmp, reducedBasisMatrix);
  }

  // add first nodes from Flux_EIM selected elems to RID
  for (int i = 1; i <= tmp.giveSize(); i++) {
    Element *elem = this->domain->giveElement(gp2elem.at(tmp.at(i)/nStressComponents)+1);
    selectedNodes.at(elem->giveNode(1)->giveNumber()) = 1;
  }

}



void
HyperReduction :: buildReducedDomain2(TimeStep *tStep, const FloatMatrix &reducedBasisMatrix) 
{

  ///testing function
  /// deactivate all elements that should be removed, deactivation of nodes is handled by numbering scheme
  IntArray nodesMultiplicity(domain->giveNumberOfDofManagers());
  nodesMultiplicity.zero();
  IntArray nonDomainNodes(domain->giveNumberOfDofManagers());
  if(!igaModeFlag) {
    this->extendReducedDomain(ridExtension, nonDomainNodes);
  } else {
    this->extendIGAReducedDomain(ridExtension, nonDomainNodes);
  }


  selectedDofs.resize(reducedBasisMatrix.giveNumberOfRows());
  selectedDofsIDMask.resize(reducedBasisMatrix.giveNumberOfRows());

  int numberOfReducedDomainDofs = 0;
  int numberOfReducedDomainNodes = 0;
  int numberOfReducedDomainElements = 0;

  int limit = max(selectedNodes.giveSize(), selectedElements.giveSize());


  for( int i = 1; i <= limit; i++) {
    if(i <= selectedNodes.giveSize()) {
      if(selectedNodes.at(i) == 1) {
	numberOfReducedDomainNodes++;
	IntArray dofIDArray;
	domain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
	for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	  if(!domain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	    numberOfReducedDomainDofs++;
	    selectedDofs.at(numberOfReducedDomainDofs) = (i-1)*dofIDArray.giveSize() + k;
	    selectedDofsIDMask.at(numberOfReducedDomainDofs) = dofIDArray.at(k);
	  }
	}
      }
    }
    
    if(i <= selectedElements.giveSize()) {
      if(selectedElements.at(i) == 1) {
	numberOfReducedDomainElements++;
	if(igaModeFlag) {
	  IGAElement *element = dynamic_cast<IGAElement*> (domain->giveElement(i));
	  if(element) {	  	  
	    for (int k = 1; k <= element->giveNumberOfIntegrationRules(); k++) {
	      IntegrationRule *iRule = element->giveIntegrationRule(k-1);
	      IGAIntegrationElement *ee = dynamic_cast< IGAIntegrationElement*> (iRule);
	      IntArray selectedIGAElementsElementI = selectedIGAElements->find(i)->second;
	      if(selectedIGAElementsElementI.at(k) != 1) {
		ee -> deActivateYourself();
	      }

	}     
	  }
	}
      
      } else {      
	this->domain->giveElement(i)->deActivateYourself();
      }
    }
  }
  

  selectedDofs.resizeWithValues(numberOfReducedDomainDofs);
  selectedDofsIDMask.resizeWithValues(numberOfReducedDomainDofs);

  
  IntArray allColumns(reducedBasisMatrix.giveNumberOfColumns());
  for(int i = 1; i <= reducedBasisMatrix.giveNumberOfColumns(); i++) {
    allColumns.at(i) = i;
  }
  hyperReducedBasisMatrix.resize(selectedDofs.giveSize(),allColumns.giveSize());
  hyperReducedBasisMatrix.beSubMatrixOf(reducedBasisMatrix,selectedDofs,allColumns);
    

  // interface nodes are : 
  // - selected nodes connected only to selected elements
  // - not Dirichlet nodes
  interfaceDofs.resize(numberOfReducedDomainDofs);
  interfaceDofs.zero();

  IntArray interfaceNodes;
  interfaceNodes.resize(numberOfReducedDomainDofs);
  interfaceNodes.zero();


  int nIDof = 0;
  int nBC = 0;
  int index = 0;    
  int nInterfaceNodes = 0;
  for(int i = 1; i <= domain->giveNumberOfDofManagers(); i++) {
    if(selectedNodes.at(domain->giveDofManager(i)->giveNumber()) == 1) {
      index++;
      IntArray dofIDArray;
      if(nonDomainNodes.at(domain->giveDofManager(i)->giveNumber()) == 1) {
	nInterfaceNodes++;
	interfaceNodes.at(nInterfaceNodes) = domain->giveDofManager(i)->giveNumber();
	domain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
	for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	  if(!domain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	    nIDof++;
	    interfaceDofs.at(nIDof) = (index-1)*dofIDArray.giveSize() + k - nBC;	
	  } else {
	    nBC++;
	  }
	}
      } else {
	domain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
	for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	  if(domain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)){
	    nBC++;	  
	  }	
	}
      }
    }
  }

  interfaceDofs.resizeWithValues(nIDof);
  interfaceNodes.resizeWithValues(nInterfaceNodes);
  //create set of interface nodes


  int nSets = domain->giveNumberOfSets() + 1;
  Set  *interfaceNodeSet = new Set(nSets, domain); //classFactory.createSet(name.c_str(), num, this)
  interfaceNodeSet->setNodeList(interfaceNodes);
  domain->resizeSets(nSets);
  domain->setSet(nSets, interfaceNodeSet);  


}




void
HyperReduction :: buildReducedDomain(TimeStep *tStep, const FloatMatrix &reducedBasisMatrix) 
{

 
  IntArray nodesMultiplicity(domain->giveNumberOfDofManagers());
  nodesMultiplicity.zero();
  IntArray nonDomainNodesTest;
  if(!igaModeFlag) {
    this->extendReducedDomain(ridExtension, nonDomainNodesTest);
  } else {
    this->extendIGAReducedDomain(ridExtension, nonDomainNodesTest);
  }


  selectedDofs.resize(reducedBasisMatrix.giveNumberOfRows());
  int numberOfReducedDomainDofs = 0;
  int numberOfReducedDomainNodes = 0;
  int numberOfReducedDomainElements = 0;

  int limit = max(selectedNodes.giveSize(), selectedElements.giveSize());



  reducedDomain = domain->Clone();

 

  for( int i = 1; i <= limit; i++) {
    if(i <= selectedNodes.giveSize()) {
      if(selectedNodes.at(i) == 1) {
	numberOfReducedDomainNodes++;
	reducedDomain->setDofManager(numberOfReducedDomainNodes, domain->giveDofManager(i));
	IntArray dofIDArray;
	domain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
	for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	  if(!reducedDomain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	    numberOfReducedDomainDofs++;
	    selectedDofs.at(numberOfReducedDomainDofs) = (i-1)*dofIDArray.giveSize() + k;
	  }
	}
      }
    }
    
    if(i <= selectedElements.giveSize()) {
      if(selectedElements.at(i) == 1) {
	numberOfReducedDomainElements++;
	reducedDomain->setElement(numberOfReducedDomainElements, domain->giveElement(i));
      }      
    }
  }
  

  selectedDofs.resizeWithValues(numberOfReducedDomainDofs);

  reducedDomain->resizeDofManagers(numberOfReducedDomainNodes);
  reducedDomain->resizeElements(numberOfReducedDomainElements);
  
  
  IntArray allColumns(reducedBasisMatrix.giveNumberOfColumns());
  for(int i = 1; i <= reducedBasisMatrix.giveNumberOfColumns(); i++) {
    allColumns.at(i) = i;
  }
  hyperReducedBasisMatrix.resize(selectedDofs.giveSize(),allColumns.giveSize());
  hyperReducedBasisMatrix.beSubMatrixOf(reducedBasisMatrix,selectedDofs,allColumns);
  
  //@todo check this
  //  reducedDomain = domain->createReducedDomain(numberOfReducedDomainNodes,selectedNodes,numberOfReducedDomainElements,selectedElements);

  

  // testing nodes are : 
  // - selected nodes connected only to selected elements
  // - not Dirichlet nodes


  IntArray nonDomainNodes(domain->giveNumberOfDofManagers());
  // this should be moved to extendRID perharps?
  for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
    if(selectedElements.at(i) == 0) {
      for (int j = 1; j <= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	nonDomainNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber()) = 1;
      }
    }    
  }

  interfaceDofs.resize(numberOfReducedDomainDofs);
  interfaceDofs.zero();



  int nIDof = 0;
  int nBC = 0;
  
  for(int i = 1; i <= reducedDomain->giveNumberOfDofManagers(); i++) {
    IntArray dofIDArray;
    if(nonDomainNodes.at(reducedDomain->giveDofManager(i)->giveGlobalNumber()) == 1) {
      reducedDomain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
      for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	if(!reducedDomain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	  nIDof++;
	  interfaceDofs.at(nIDof) = (i-1)*dofIDArray.giveSize() + k - nBC;	
	} else {
	  nBC++;
	}
      }
    } else {
      reducedDomain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
      for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	if(reducedDomain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)){
	  nBC++;	  
	}	
      }
    }
  }

  interfaceDofs.resizeWithValues(nIDof);


}


void 
HyperReduction :: extendReducedDomain(int times, IntArray &nonDomainNodes)
{

  nonDomainNodes.resize(domain->giveNumberOfDofManagers());
  nonDomainNodes.zero();

  IntArray element2ReduceNumbers;
  /// one set containing elements to be always reduced
  if(ridElement2ReduceSetNumbers.giveSize()) {
    int setId = ridElement2ReduceSetNumbers.at(1);
    element2ReduceNumbers = this->domain->giveSet(setId)->giveElementList();
  }
  
  // add connected elements to the selected ones
  for (int count = 1; count <= times; count++) {
    
    // select the elements attached to selected nodes
    for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
      for (int j = 1; j<= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	if (selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber())) {
	  selectedElements.at(i) = 1; 
	  break;
	}
      }
    }
    
    /// reduced elements which are in the element2reduceset
    for(int i = 1; i <= element2ReduceNumbers.giveSize(); i++) {
      selectedElements.at(element2ReduceNumbers.at(i)) = 0;

    }


    
    // select the nodes attached to selected elements
    for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
      if(! selectedElements.at(i))
	continue;
      for (int j = 1; j <= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber()) = 1;
      }
    }
    
  }


  for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
    if(selectedElements.at(i) == 0) {
      for (int j = 1; j <= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	nonDomainNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber()) = 1;
      }
    }    
  }



  
}
 

void 
HyperReduction :: extendIGAReducedDomain(int times,  IntArray &nonDomainNodes)
{



  this->selectedIGAElements = new std :: map< int, IntArray >();

  // add connected elements to the selected ones
  for (int count=0; count<times; count++) {
    if(count == times-1) {
      // this was moved here from computeReducedDomain_Nodes
      for(int i = 1; i <= ridNodeSetNumbers.giveSize(); i++) {
	int setId = ridNodeSetNumbers.at(i);
	IntArray nodesNumbers = this->domain->giveSet(setId)->giveNodeList();
	for(int j = 1; j <= nodesNumbers.giveSize(); j++) {
	  selectedNodes.at(nodesNumbers.at(j)) = 1;
	}
      }
      ///////////////////////////////////////////////////////////
    }


    IntArray selectedIGAElementsElementI;
    // select the elements attached to selected nodes
    for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
      for (int j = 1; j<= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	if (selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber())) {
	  selectedElements.at(i) = 1; 
	  ///// IGA	  
	  IGAElement *element = dynamic_cast<IGAElement*> (domain->giveElement(i));
	  if(element) {	  	  
	    selectedIGAElementsElementI.resize(element->giveNumberOfIntegrationRules());
	    selectedIGAElementsElementI.zero();
	    IntArray mask;
	    for (int k = 1; k <= element->giveNumberOfIntegrationRules(); k++) {
	      IntegrationRule *iRule = element->giveIntegrationRule(k-1);
	      IGAIntegrationElement *ee = dynamic_cast< IGAIntegrationElement*> (iRule);
	      element->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
	      for(int m = 1; m <= mask.giveSize(); m++) {
		int number = domain->giveElement(i)->giveNode(mask.at(m))->giveNumber();
		if(selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(m))->giveNumber()) == 1){
		   //selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber()))
		  selectedIGAElementsElementI.at(k) = 1;
		  break;
		}
	      }
	    }
	  }
	  (*selectedIGAElements)[i] = selectedIGAElementsElementI;
	  /////////////////////////////////
	  break;
	}
      }
    }
    
    
    // select the nodes attached to selected IGA elements
    for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
      if(!selectedElements.at(i))
	continue;
      
      IGAElement *element = dynamic_cast<IGAElement*> (domain->giveElement(i));
      if(element) {	  	  
	IntArray selectedIGAElementsElementI;
	selectedIGAElementsElementI = selectedIGAElements->find(i)->second;
	IntArray mask;
	for (int k = 1; k <= element->giveNumberOfIntegrationRules(); k++) {
	  if(selectedIGAElementsElementI.at(k) == 1) {
	    IntegrationRule *iRule = element->giveIntegrationRule(k-1);
	    IGAIntegrationElement *ee = dynamic_cast< IGAIntegrationElement*> (iRule);
	    element->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
	    for (int j = 1; j <= mask.giveSize(); j++) {
	      if(count == times-1) { // last iteration of extension, take these nodes as nondomain nodes, so they will be in the set of interface nodes
	
		if( selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) == 1) {
		} else {
		  nonDomainNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) = 1;
		  selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) = 1;
		}		  
	      }	else {
		selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) = 1;
	      }
	    }
	  }
	}
      }
    }
  }
}



  /*
void 
HyperReduction :: extendIGAReducedDomain2(int times,  IntArray &nonDomainNodes)
{



  this->selectedIGAElements = new std :: map< int, IntArray >();
  // add connected elements to the selected ones
  
  // first add nodes given by node set
  for(int i = 1; i <= ridNodeSetNumbers.giveSize(); i++) {
    int setId = ridNodeSetNumbers.at(i);
    IntArray nodesNumbers = this->domain->giveSet(setId)->giveNodeList();
    for(int j = 1; j <= nodesNumbers.giveSize(); j++) {
      selectedNodes.at(nodesNumbers.at(j)) = 1;
    }
  }
  
  IntArray selectedIGAElementsElementI;
  // select the elements attached to selected nodes
  for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
    for (int j = 1; j<= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
      if (selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber())) {
	selectedElements.at(i) = 1; 
	///// IGA	  
	IGAElement *element = dynamic_cast<IGAElement*> (domain->giveElement(i));
	if(element) {	  	  
	  selectedIGAElementsElementI.resize(element->giveNumberOfIntegrationRules());
	  selectedIGAElementsElementI.zero();
	  IntArray mask;
	  for (int k = 1; k <= element->giveNumberOfIntegrationRules(); k++) {
	    IntegrationRule *iRule = element->giveIntegrationRule(k-1);
	    IGAIntegrationElement *ee = dynamic_cast< IGAIntegrationElement*> (iRule);
	    element->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
	    for(int m = 1; m <= mask.giveSize(); m++) {
	      int number = domain->giveElement(i)->giveNode(mask.at(m))->giveNumber();
	      if(selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(m))->giveNumber()) == 1){
		//selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber()))
		selectedIGAElementsElementI.at(k) = 1;
		break;
	      }
	      }
	  }
	}
	(*selectedIGAElements)[i] = selectedIGAElementsElementI;
	/////////////////////////////////
	break;
      }
    }
  }
    
    
  // select the nodes attached to selected IGA elements
  for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
    if(!selectedElements.at(i))
      continue;
    IGAElement *element = dynamic_cast<IGAElement*> (domain->giveElement(i));
    if(element) {	  	  
      IntArray selectedIGAElementsElementI;
      selectedIGAElementsElementI = selectedIGAElements->find(i)->second;
      IntArray mask;
      for (int k = 1; k <= element->giveNumberOfIntegrationRules(); k++) {
	if(selectedIGAElementsElementI.at(k) == 1) {
	  IntegrationRule *iRule = element->giveIntegrationRule(k-1);
	  IGAIntegrationElement *ee = dynamic_cast< IGAIntegrationElement*> (iRule);
	  element->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
	  for (int j = 1; j <= mask.giveSize(); j++) {
	    if(count == times-1) { // last iteration of extension, take these nodes as nondomain nodes, so they will be in the set of interface nodes
	      if( selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) == 1) {
		} else {
		  nonDomainNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) = 1;
		  selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) = 1;
		}		  
	      }	else {
		selectedNodes.at(domain->giveElement(i)->giveNode(mask.at(j))->giveNumber()) = 1;
	      }
	    }
	  }
	}
      }
    }
  }
}
  */

bool 
HyperReduction :: EIM(IntArray &RID_ddl, const FloatMatrix &reducedBasisMatrix) 
{
  int rk_max; 
  FloatArray Resid, tmp;
  FloatMatrix V;
  // V = A(:,1)/A(imax,1); RID_ddl = [RID_ddl,imax];  
  Resid.beColumnOf(reducedBasisMatrix,1);
  V.resize(Resid.giveSize(), 1);
  V.setColumn(Resid, 1);
  //  Resid.beColumnOf(V,1);
  rk_max = Resid.giveIndexMaxAbsElem(); 
  Resid.times(1./fabs(Resid.at(rk_max)));
  // put index of maximal value to the RID_ddl
  RID_ddl.resizeWithValues(RID_ddl.giveSize()+1);
  RID_ddl.at(RID_ddl.giveSize()) = rk_max;
  
  FloatArray gamma, b;
  FloatMatrix VDDL, tmpMtrx;
  IntArray kArray;
  for (int k=2; k<= reducedBasisMatrix.giveNumberOfColumns(); k++) {
    // gamma = V(RID_ddl,:)\A(RID_ddl,k);
    kArray.resize(1);
    kArray.at(1) = k;
    tmpMtrx.beSubMatrixOf(reducedBasisMatrix, RID_ddl, kArray);
    b.beColumnOf(tmpMtrx, 1);

    IntArray lArray, rArray;
    lArray.resize(V.giveNumberOfColumns());
    rArray.resize(V.giveNumberOfColumns());
    for(int i = 1; i <= V.giveNumberOfColumns(); i ++) {
      lArray.at(i) = RID_ddl.at(i);
      rArray.at(i) = i;
    }
    VDDL.resize(V.giveNumberOfColumns(),V.giveNumberOfColumns());
    VDDL.beSubMatrixOf(V,lArray,rArray);
    // VDDL.beSubMatrixOf(V, RID_ddl.at(1),RID_ddl.at(RID_ddl.giveSize()) , 1, V.giveNumberOfColumns());
    // solve gamma = b/VDDL;
    VDDL.solveForRhs(b, gamma);    
    // R = A(:,k) - V*gamma;
    Resid.beColumnOf(reducedBasisMatrix,k);
    tmp.beProductOf(V,gamma);
    tmp.times(-1.);
    Resid.add(tmp);

    // [val,imax] = max(abs(R)); RID_ddl = [RID_ddl,imax]; V = [V, R/R(imax)];
    rk_max = Resid.giveIndexMaxAbsElem(); 
    Resid.times(1./fabs(Resid.at(rk_max)));    
    // put index of maximal value to the RID_ddl
    RID_ddl.resizeWithValues(RID_ddl.giveSize()+1);
    RID_ddl.at(RID_ddl.giveSize()) = rk_max;
    V.resizeWithData(V.giveNumberOfRows(), V.giveNumberOfColumns()+1);
    V.setColumn(Resid, V.giveNumberOfColumns());
    
  }
  
  return true;
}








} // end namespace oofem
