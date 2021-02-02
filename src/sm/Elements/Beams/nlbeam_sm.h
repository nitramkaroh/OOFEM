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

#ifndef nlbeam_sm_h
#define nlbeam_sm_h


#include "dofmanager.h"
#include "Elements/nlstructuralelement.h"


///@name Input fields for NlBeam_SM
//@{
#define _IFT_NlBeam_SM_Name "nlbeam_sm"
#define _IFT_NlBeam_SM_NIP "nip"
#define _IFT_NlBeam_SM_EA "ea"
#define _IFT_NlBeam_SM_EI "ei"
#define _IFT_NlBeam_SM_Beam_Tolerance "btol"
#define _IFT_NlBeam_SM_Beam_MaxIteration "bmaxit"

//@}

namespace oofem {


/**
 * This class implements a 2-dimensional large strain beam element
 * The shooting method is used to calculate internal forces and stiffness matrix
 * Add more description 
 */
class NlBeam_SM : public NLStructuralElement
{
protected:
    int NIP = 100;
    double pitch = 10, beamLength = 0;
    FloatArray internalForces;
    FloatArray x, u, w, phi;
    FloatMatrix jacobi;
    double beam_tol = 1.e-6, beam_maxit = 100;
    double EI, EA;
    FloatArray vM, vV, vN;
public:
    NlBeam_SM(int n, Domain *aDomain);
    virtual ~NlBeam_SM(){;}

    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    virtual int computeNumberOfDofs() { return 6; }

    virtual void  printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "NlBeam_SM"; }
    virtual const char *giveInputRecordName() const { return _IFT_NlBeam_SM_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    

protected:

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
			  TimeStep *tStep = NULL, int lowerIndx = 1, int upperIndx = ALL_STRAINS){;}
    virtual double computeLength();
    double givePitch();
    virtual MaterialMode giveMaterialMode() { return _2dBeam; }

    double computeMomentFromCurvature(double kappa);
    double computeDerMomentFromCurvature(double kappa);
    double computeCurvatureFromMoment(double M);
    void   integrateAlongBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi);
    void findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc);
    void construct_T(FloatMatrix &T, const double phia);
    void construct_Tprime(FloatMatrix &T, const double phia);
    void construct_l(FloatArray &l, double phia);
    void construct_l(FloatArray &l, double phia, double L);
      
    void construct_lprime(FloatArray &l, const double phia);
    void findLeftEndForces(const FloatArray &u, FloatArray &fab);
    
    FILE * giveOutputStream(TimeStep *tStep);


    void giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u);
    void computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

};
} // end namespace oofem
#endif // beam2d_h
