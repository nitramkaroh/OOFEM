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
 *               Copyright (C) 1993 - 2014    Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef fbarelementinterface_h
#define fbarelementinterface_h

#include "interface.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "nlstructuralelement.h"



///@name fbarelementinterface
//@{
#define _IFT_FbarElementExtensionInterface_fbarflag "fbarflag"
//@}

namespace oofem {

class GaussPoint;
class TimeStep;

/**
 * F-bar element extension interface
 * @author Martin Horak, nitramkaroh@seznam.cz
 * @note Reference: EA de Souze Neto, D. Peric, M. Dutko, DRJ Owen: Design of simple low order finite elements for large strain analysis of nearly incompressible solids
 * The method is based on replacing of deformation gradient by its modified counterpart
 * \bar{F} = (\frac{J_0}{J})^{/frac{1}{3}}F
 * Appropriate to tackle incompressibility locking in large strains
 * Elements using F-bar interface just have to be inherited from this interface : Note that this approach is not appropriate for triangles and tetras
 * 
 */

class FbarElementExtensionInterface : public Interface
{



protected:
  /**
   * Flag indicating whether F-bar formulation is used
   */
  int FbarFlag;
  int initialStepFlag;
  FloatArray u;
		
  
public:
    /**
     * Constructor. Creates element interface belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    FbarElementExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~FbarElementExtensionInterface() { } 
  
    IRResultType initializeFrom(InputRecord *ir);


    /**
     * Returns FbarFlag
     */
    int giveFbarFlag(){return FbarFlag;}
    int giveInitialStepFlag(){return initialStepFlag;}
    void setIntialStepFlag(int isf){initialStepFlag = isf;}

    FloatArray &giveU(){return u;}
    void setU(FloatArray &disp){u = disp;}



    virtual void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem, ValueModeType vmt = VM_Total);


    /**
     * Computes the modified deformation gradient Fbar in Voigt form at integration point ip and at time
     * step tStep. Computes the displacement gradient and adds an identitiy tensor.
     *
     * @param answer Deformation gradient vector
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    virtual double computeFbarDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem, ValueModeType vmt = VM_Total);

    /**
     * Computes the stiffness matrix corrected due to F-bar formulation.
     * The response is evaluated using @f$ \int B_{\mathrm{H}}^{\mathrm{T}} D (B_{\mathrm{H}})_0 \;\mathrm{d}v @f$, where
     * @f$ B_{\mathrm{H}} @f$ is the B-matrix which produces the displacement gradient vector @f$ H_{\mathrm{V}} @f$ when multiplied with
     * the solution vector a.
     * @f $ (B_{\mathrm{H}})_0 @f$ is the B-matrix which produces the displacement gradient vector in the element parametric centroid
     * @param answer Computed stiffness matrix correction.
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    virtual void computeFbarStiffnessMatrix(FloatMatrix &answer1, MatResponseMode rMode, TimeStep *tStep, NLStructuralElement *e );
    

    /**
     * Computes two correction terms due to the F-bar formulation
     * @f$ \bar{A}^1_{ijmn} = (\frac{J_0}{J})^{1/3}(A_{ijmn} - \frac{1}{3}A_ijkl F_{nm}^{-1}F_{kl}) @f
     * @f$ \bar{A}^2_{ijmn} = \frac{1}{3}(\frac{J_0}{J})^{1/3}(A_ijkl (F_{nm})^{-1}_0F_{kl}) @f
     * @param answer Computed tranformed first elasticty stiffness.
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    virtual double computeFbarStiffnessCorrections(FloatMatrix &answer1, FloatMatrix &answer2, const FloatMatrix &dPdF, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem);

    virtual void computeSpatialFbarStiffnessCorrections(FloatMatrix &answer, FloatMatrix &a, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem);

 
};
} 
#endif //fbarelementinterface_h
