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

#ifndef nlbeam_internalcontact_h
#define nlbeam_internalcontact_h

#include "Elements/Beams/nlbeam_sm.h"


///@name Input fields for NlBeam_SM2
//@{
#define _IFT_NlBeamInternalContact_Name "nlbeam_internalcontact"

//@}

namespace oofem {
 enum contactModeType {N_cmode, AA_cmode, AB_cmode, AC_cmode, BA_cmode, BB_cmode, BC_cmode, CA_cmode, CB_cmode};
  enum processType {Unknown_proc, Roll_proc, Stick_proc, Slid_proc};


/**
 * This class implements a 2-dimensional large strain beam element with internal contact
 * The shooting method is used to calculate internal forces and stiffness matrix
 * Add more description 
 */
class NlBeamInternalContact : public NlBeam_SM
{
protected:
  // coefficient of friction
  double friction =  0.3;
  double  DX  = 0.11;
  // tolerance for iterations at the beam level (negligible rotation)
  double TOL_BEAM =  1.e-8;
  // relative tolerance for beam length differences to detect sliding
  double TOL_LENGTH =  1.e-8;
  // maximum number of iterations at the beam level
  int MAXIT_BEAM  = 20;       
  // maximum number of iterations of the assumed contact mode
  int MAXIT_CONTACT_MODE = 5;
  // tolerance for iterations of the contact time
  double TOL_CONTACT_TIME =  1.e-8;
  // maximum number of iterations of the contact time
  int MAXIT_CONTACT_TIME =  10;

  double leftSegmentLength = 60.;//40.;
  double rightSegmentLength = 0.001;//15.;
  double leftActiveSegmentLength = 49.999;//35.;
  double rightActiveSegmentLength = 0.001;//15.;
  double trialLeftActiveSegmentLength, trialRightActiveSegmentLength;
  
 
#define PI 3.1415926
  
  double Jacobi[3][3], Jacobi44[4][4], Kblock[3][3];
  double alpha = 0;
  contactModeType contactMode = AC_cmode;
  contactModeType trialContactMode;
  processType Process = Unknown_proc;
  processType trialProcess = Unknown_proc;
  bool stiffEvalMode = false;
       
public:
    NlBeamInternalContact(int n, Domain *aDomain);
    virtual ~NlBeamInternalContact(){;}

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual const char *giveClassName() const { return "NlBeamInternalContact"; }
    virtual const char *giveInputRecordName() const { return _IFT_NlBeamInternalContact_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

protected:

    
  double signum(double x);
  double L2norm(double x, double y, double z);
  void solve(double A[4][4], double b[4], double x[4], int n);
  void solve3(double A[3][3], double b[3], double x[3]);
  void invert(double A[3][3], double Ainv[3][3]);
  void invert4(double A[4][4], double Ainv[4][4]);
  double computeCurvatureFromMoment(double M);
  double computeDerMomentFromCurvature(double kappa);
  double evalContactLoadingFunction(double Nc, double Qc, contactModeType cmode);
  double evalDerContactLoadingFunction(double Nc, double dN, double dQ, contactModeType cmode);
  bool integrateAlongSegment(double fab[3], double deltaPhi, double Lb, double ub[3], double jac_b[3][4], double* Lc, double uc[3], double jac_c[3][4], bool inflection_only, FILE* outfile);
  void transform_ub2ua(double ub[3], double ua[3]);
  void transform_fab2fba(double ub[3], double fab[3], double fba[3]);
  void transform_fba2fab(double ub[3], double fba[3], double fab[3]);
  bool checkLeftSegmentLengthAdmissibility(double Lc);
  bool checkRightSegmentLengthAdmissibility(double Lc);
  bool checkRotationAdmissibility(double deltaPhi, contactModeType cmode);
  bool findLeftEndForcesLocal_Tip_SoS(double ub_target[3], double fab[3], double* Lac, double Lb, contactModeType cmode, int process, double* Nca, double* Qca, double* deltaPhi, bool printflag);
  double shiftToIntervalFromMinusPiToPi(double x);
  contactModeType suggestedMode(double deltaPhi, bool seglength_ok, contactModeType current_cmode);
  contactModeType findLeftEndForcesLocal_Tip(double ub_target[3], double fab[3], contactModeType cmode, bool printflag);
  bool findLeftEndForcesLocal_Smooth_Rolling(double ub_target[3], double fab[3], double deltaPhi, double* Lac, double* Lbc, double* Nc, double* Qc, bool printflag);
  bool findLeftEndForcesLocal_Smooth_Sliding(double ub_target[3], double fab[3], double deltaPhi, double* Lac, double* Lbc, contactModeType cmode, bool printflag);
  contactModeType findLeftEndForcesLocal_Smooth(double ub_target[3], double fab[3], contactModeType cmode, bool printflag);
  contactModeType predictContactMode(double ub[3]);
  double findTipContactTime(double u_prev[3], double du[3], double segLength);
  contactModeType predictContactMode(double ub[3], double ub_prev[3]);
  bool findLeftEndForcesLocal(double ub[3], double ub_prev[3], double fab[3], bool printflag);
  void construct_T(double T[3][3], double phia);
  void construct_Tprime(double T[3][3], double phia);
  void construct_l(double l[3], double phia);
  void construct_l_IC(double l[3], double phia, double L);
  void construct_lprime(double l[3], double phia);
  bool findLeftEndForces(double u[6], double u_prev[6], double fab[3]);
  bool findEndForces(double u[6], double u_prev[6], double f[6]);
  bool evalStiffnessMatrix(double u[6], double u_prev[6], double fab[3], double K[6][6]);
  void integrateAlongSegmentAndPlot(double fab[3], double Lb, double segmentLength, double u0[2], double T[2][2], FILE* outfile);
  void plotSegment(double fab[3], double ub[3], bool isLeftSegment, FILE* outfile);
  void plotResponse(int nstage, int nstep[10], double ustep[3][10], int iplot[10]);
};
} // end namespace oofem
#endif // nlbem_internalcontact_h
