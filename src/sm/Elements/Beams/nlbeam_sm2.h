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

#ifndef nlbeam_sm2_h
#define nlbeam_sm2_h


#include "dofmanager.h"
#include "Elements/nlstructuralelement.h"
#include "scalarfunction.h"


///@name Input fields for NlBeam_SM2
//@{
#define _IFT_NlBeam_SM2_Name "nlbeam_sm2"
#define _IFT_NlBeam_SM2_NIP "nip"
#define _IFT_NlBeam_SM2_EA "ea"
#define _IFT_NlBeam_SM2_EI "ei"
#define _IFT_NlBeam_SM2_initialCurvature "initialcurvature"
#define _IFT_NlBeam_SM2_RADIUS "radius"
#define _IFT_NlBeam_SM2_DEPTH "depth"
#define _IFT_NlBeam_SM2_Material "materialtype"
#define _IFT_NlBeam_SM2_Beam_Tolerance "btol"
#define _IFT_NlBeam_SM2_Beam_MaxIteration "bmaxit"
#define _IFT_NlBeam_SM2_Section_Tolerance "stol"
#define _IFT_NlBeam_SM2_Section_MaxIteration "smaxit"

#define _IFT_NlBeam_SM2_fu0 "fu0"
#define _IFT_NlBeam_SM2_fw0 "fw0"
#define _IFT_NlBeam_SM2_fphi0 "fphi0"
#define _IFT_NlBeam_SM2_fkappa0 "fkappa0"
#define _IFT_NlBeam_SM2_tangentPoint "tp"

#define _IFT_NlBeam_SM2_alCoord "alcoord"
#define _IFT_NlBeam_SM2_pu0 "pu0"
#define _IFT_NlBeam_SM2_pw0 "pw0"
#define _IFT_NlBeam_SM2_pphi0 "pphi0"
#define _IFT_NlBeam_SM2_pkappa0 "pkappa0"
#define _IFT_NlBeam_SM2_beta "beta"
#define _IFT_NlBeam_SM2_curvedbeamLength "curvedbeamlength"


//@}

namespace oofem {


/**
 * This class implements a 2-dimensional large strain beam element
 * The shooting method is used to calculate internal forces and stiffness matrix
 * Add more description 
 */
class NlBeam_SM2 : public NLStructuralElement
{
protected:
    int NIP = 100;
    double pitch = 10, beamLength = 0;
    FloatArray internalForces;
    FloatArray x, u, w, phi;
    FloatMatrix jacobi;
    double beam_tol = 1.e-6, beam_maxit = 100;
    double section_tol = 1.e-6,section_maxit = 20;
    double EI, EA;
    double RADIUS, DEPTH;
    double curvedbeamLength;
    double choordLength;
    int materialtype = 1;
    FloatArray vM, vV, vN;
    double cosBeta, sinBeta, cosAlpha, sinAlpha;
    /// Initial stress field curved configuration
    ScalarFunction fu_0, fw_0, fphi_0, fkappa_0;
    FloatArray pu_0, pw_0, pphi_0, pkappa_0;
    FloatArray u_0, w_0, phi_0, kappa_0, phi_12;
    double beta = 0;
    FloatArray tangentPoint;
    FloatArray arcLengthCoordinate;

    
public:
    NlBeam_SM2(int n, Domain *aDomain);
    virtual ~NlBeam_SM2(){;}

    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    virtual int computeNumberOfDofs() { return 6; }

    virtual void  printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "NlBeam_SM2"; }
    virtual const char *giveInputRecordName() const { return _IFT_NlBeam_SM2_Name; }
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
    void integrateAlongBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi);
    void integrateAlongStraightBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi);
    void integrateAlongCurvedBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi);
    void findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc);
    void construct_T(FloatMatrix &T, const double phia);
    void construct_Tprime(FloatMatrix &T, const double phia);
    void construct_l(FloatArray &l, double phia);
    void construct_l(FloatArray &l, double phia, double L);
    void construct_lprime(FloatArray &l, const double phia);
  
    void findLeftEndForces(const FloatArray &u, FloatArray &fab);

    void  printOutputAt_StraightBeam(FILE *file, TimeStep *tStep);
    void  printOutputAt_CurvedBeam(FILE *file, TimeStep *tStep);

    double eval_kappa0(double x);
    double eval_phi0(double x);
    double eval_u0(double x);
    double eval_w0(double x);
    
    double eval_c1(FloatArray &u);
    double eval_c2(FloatArray &u);
    double eval_c1_StraightBeam(FloatArray &u);
    double eval_c2_StraightBeam(FloatArray &u);
    double eval_c1_CurvedBeam(FloatArray &u);
    double eval_c2_CurvedBeam(FloatArray &u);
    double computeDeltaCurvatureFromInternalForces(double M, double N, double curvature);
    double computeDerCurvatureMoment(double M, double N, double curvature);
    double computeDerCurvatureNormalForce(double M, double N, double curvature);
    double computeCenterlineStrainFromInternalForces(double M, double N, double curvature);
    double computeDerStrainMoment(double M, double N, double curvature);
    double computeDerStrainNormalForce(double M, double N);

    FILE * giveOutputStream(TimeStep *tStep);


    void computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);    
    void giveInternalForcesVector_from_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u);


};
} // end namespace oofem
#endif // beam2d_h
