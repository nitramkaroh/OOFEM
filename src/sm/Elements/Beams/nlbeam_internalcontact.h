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
#define _IFT_NlBeamInternalContact_LeftSegmentLength "lsl"
#define _IFT_NlBeamInternalContact_RightSegmentLength "rsl"
#define _IFT_NlBeamInternalContact_LeftActiveSegmentLength "lasl"
#define _IFT_NlBeamInternalContact_RightActiveSegmentLength "rasl"
#define _IFT_NlBeamInternalContact_ContactMode "cmode"
#define _IFT_NlBeamInternalContact_Friction "friction"
#define _IFT_NlBeamInternalContact_dx "dx"

//@}

namespace oofem {
 enum ContactModeType {CMT_N, CMT_AA, CMT_AB, CMT_AC, CMT_BA, CMT_BB, CMT_BC, CMT_CA, CMT_CB};
 enum ProcessType {PT_Unknown, PT_Roll, PT_Stick, PT_Slide};

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

  double leftSegmentLength = 40.;
  double rightSegmentLength = 0.001;//15.;
  double leftActiveSegmentLength = 49.999;//35.;
  double rightActiveSegmentLength = 0.001;//15.;
  double trialLeftActiveSegmentLength, trialRightActiveSegmentLength;
  FloatArray internalForces;
 
  FloatMatrix Jacobi;
  FloatMatrix Jacobi44;
  FloatMatrix Kblock;
  double alpha = 0;
  ContactModeType contactMode = CMT_AC;
  ContactModeType trialContactMode;
  ProcessType Process = PT_Unknown;
  ProcessType trialProcess = PT_Unknown;
  bool stiffEvalMode = false;
  
       
public:
    NlBeamInternalContact(int n, Domain *aDomain);
    virtual ~NlBeamInternalContact(){;}
    
    /*
      Evaluate the tangent stiffness matrix based on given displacements at both ends
      and given end forces at the left end (must be provided).
      It is assumed that findLeftEndForces has been run first and that the 3x3 Jacobi matrix
      has been stored in Jacobi[3][3] or the 4x4 Jacobi matrix has been stored in Jacobi44[4][4].
      
      Input variables:
      u[6] ... end displacements and rotations at the end of the current step, in global coordinates
      u_prev[6] ... end displacements and rotations at the end of the previous step, in global coordinates
      fab[3] ...  a guess (or converged value) of the left-end forces and moment, in local coordinates
      
      Output variables:
      fab[3] ...  left-end forces and moment (values of initial guess are rewritten)
      K[6][6] ... element stiffness matrix in global coordinates
      
      Return value:
      boolean, indicating whether everything has run smoothly
    */
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual const char *giveClassName() const { return "NlBeamInternalContact"; }
    virtual const char *giveInputRecordName() const { return _IFT_NlBeamInternalContact_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    void giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );
  Interface *giveInterface(InterfaceType it) override;
     Element_Geometry_Type giveGeometryType() const override { return EGT_Composite; }


    
protected:

    
  double L2norm(double x, double y, double z);
  double computeCurvatureFromMoment(double M);
  double computeDerMomentFromCurvature(double kappa);

  /*
    Evaluation of the contact loading function fc from the given normal force and shear force 
    at the contact point. Negative values of fc indicate sticking and zero value indicates sliding.
    The shear force is actually transmitted in the direction normal to the contact surface and
    the normal force in the direction tangential to the contact surface. Whether positive shear
    force means tension (inadmissible) or compression in the contact surface depends on which
    segment is above and which one is below. Positive shear force means that the left segment
    acts on the right segment upwards.
    
    Input:
    Nc ... normal force at contact point
    Qc ... shear force at contact point
    cmode ... contact mode
    
    Return value:
    fc ... contact loading function
  */   
  double evalContactLoadingFunction(double Nc, double Qc, ContactModeType cmode);
  /*
    Evaluation of the differential (linearized increment) of the contact loading function fc.
    
    Input:
    Nc ... normal force at contact point
    dN ... increment of normal force at contact point
    dQ ... increment of shear force at contact point
    cmode ... contact mode
    
    Return value:
    dfc ... linearized increment of the contact loading function
  */   
  double evalDerContactLoadingFunction(double Nc, double dN, double dQ, ContactModeType cmode);
  // ========================================================================
// KEY ALGORITHM - INTEGRATION ALONG A BEAM SEGMENT
// ========================================================================

/*
  Numerical integration of differential equations describing a beam segment.
  Left-end displacements set to zero, left-end forces specified as input.
  Computes displacements at the end of the segment of given length Lb.
  If an inflection point c is detected, computes its distance from the left end, Lc, and the displacements at c.
  Also computes two Jacobi matrices that contain derivatives of displacements at b or c with respect to the left-end forces
  and with respect to the spatial coordinate. 
  If the user specifies it, integration stops at the first inflection point and the results that refer to the segment end are then meaningless.
  Everything is done in local coordinates aligned with the left beam end.

  Input variables:
  * fab ... horizontal force, vertical force and moment at the left end
  * deltaPhi ... prescribed displacement jump at the inflexion point (useful for inverted modes, otherwise zero) 
  * Lb ... total length of the segment
  * inflection_only ... if this is set to "true", the integration stops as soon as an inflection point is found

  Output variables:
  * ub ... horizontal displacement, vertical displacement and rotation at the end of the segment
  * jac_b ... derivatives of ub with respect to fab and x (evaluated at x=Lb)
  * Lc ... coordinate of inflection point
  * uc ... horizontal displacement, vertical displacement and rotation at the inflection point
  * jac_c ... derivatives of uc with respect to fab and x (evaluated at x=Lc)

  Return value:
  * true or false, indicating whether an inflection point has been detected (if not, the output Lc, uc and jac_c is meaningless)
*/
  bool integrateAlongSegment(FloatArray &fab, double deltaPhi, double Lb, FloatArray &ub, FloatMatrix &jac_b, double &Lc, FloatArray &uc, FloatMatrix &jac_c, bool inflection_only);
  /*
  This method computes the right-end forces and moment from the left-end forces and moment
  and from the relative displacements and rotation of the right end, based on equilibrium.
  
  Input: 
  ub ... right-end relative displacements and rotation
  fab ... left-end forces and moment
  
  Output: 
  fba ... right-end forces and moment
*/ 

  void transform_ub2ua(const FloatArray &ub, FloatArray &ua);
  /*
  This method computes the right-end forces and moment from the left-end forces and moment
  and from the relative displacements and rotation of the right end, based on equilibrium.
  
  Input: 
  ub ... right-end relative displacements and rotation
  fab ... left-end forces and moment
  
  Output: 
  fba ... right-end forces and moment
*/ 
  void transform_fab2fba(const FloatArray &ub, const FloatArray &fba, FloatArray &answer);
  /*
  This method computes the right-end forces and moment from the left-end forces and moment
  and from the relative displacements and rotation of the right end, based on equilibrium.

  Input: 
  ub ... right-end relative displacements and rotation
  fba ... right-end forces and moment
  
  Output: 
  fab ... left-end forces and moment
*/ 
  void transform_fba2fab(const FloatArray &ub, const FloatArray &fba, FloatArray &answer);
  /*
   */
  
  bool checkLeftSegmentLengthAdmissibility(double Lc);
  bool checkRightSegmentLengthAdmissibility(double Lc);
  /*
  For the given contact mode (which must be one of the tip modes),
  this method checks whether the given rotation jump at the contact point would be admissible.

  Input:
  deltaPhi ... rotation jump
  cmode ... contact mode (one of the tip modes)
  
  Return value:
  'true' if the rotation jump is admissible
*/
  bool checkRotationAdmissibility(double deltaPhi, ContactModeType cmode);
  /*
 Iterative search for the left-end forces and variables that characterize the contact behavior.
 Everything is done in local coordinates aligned with the left beam end.
 The tip of the right segment is assumed to be in contact with some point on the left segment.
 For modes AC and BC, this assumption is satisfied naturally.
 For modes CA and CB, it is necessary to swap left and right, which is done before the present method is invoked.
 It is specified whether "sticking" is assumed, in which case the position of the contact point is prescribed 
 and does not change, or "sliding" is assumed, in which case this position needs to be determined.
 Note that this function does NOT check admissibility of the solution. It only makes sure that
 three conditions are satisfied. Two of them are compatibility conditions (in the deformed
 configuration, the position of the contact point on the left segment coincides with the
 position of the tip of the right segment). The third condition is "zero moment at the (fixed)
 contact point" for the sticking mode and "zero contact loading function at the contact point"
 for the sliding mode. The loading function is computed from the normal and shear forces and
 its zero value indicates the critical combination that leads to frictional sliding.

 Input variables: 
 * ub_target ... displacements and rotation of the right end wrt to the fixed left end that should be achieved
 * fab ... initial guess of the left-end forces 
 * Lac ... initial active length of the left segment (can change if sliding occurs)
 * Lb ... total length of the right segment (is fixed, tip contact assumed)
 * cmode ... contact mode (one of the tip modes)
 * process ... type of process (0=sticking, 1=sliding)
 * printflag ... flag indicating whether detailed results should be printed

 Output variables: 
 * fab ... left-end forces at the end of the step
 * Lac ... final active length of the left segment (different from initial value if sliding occured)
 * Nca ... normal force in the left segment at the contact point (= inflection point)
 * Qca ... shear force in the left segment at the contact point
 * deltaPhi ... relative rotation of the right segment wrt to the left segment at the contact point
 
 Return value:
 * true or false, indicates whether a converged solution has been found (admissibility is not checked here)
*/
  bool findLeftEndForcesLocal_Tip_SoS(const FloatArray &ub_target, FloatArray &fab, double &Lac, double Lb, ContactModeType cmode, int process, double &Nca, double &Qca, double &deltaPhi);
  /* 
   Auxiliary function: the argument 'x' is shifted by an integer multiple of 2.*PI to get a number between -PI and PI.
*/
  double shiftToIntervalFromMinusPiToPi(double x);
  /*
  The most likely contact mode at the end of the step is determined based on the information
  on the results obtained with the assumption that the mode is equal to a given tip mode.
  The information provided is the rotation jump at the contact point and the result of
  segment length admissibility evaluation.

  Input:
  deltaPhi ... rotation jump at the contact point
  seglength_ok ... 'true' if the active length of both segments is between zero and the total segment length
  current_cmode ... the mode for which the solution has been computed (must be one of the tip modes)

  Return value:
  The contact mode that is expected to occur at the end of the step.
*/
  ContactModeType suggestedMode(double deltaPhi, bool seglength_ok, ContactModeType current_cmode);
  /*
 Evaluation of the left-end forces and moment for given relative displacements and rotation
 of the right end, assuming tip contact.
 Everything is done in local coordinates aligned with the left beam end.
 The tip of one segment is assumed to be in contact with some point on the other segment.
 For mode AC, the tip of the right segment touches the left segment from "above", for mode BC from "below".
 For mode CB, the tip of the left segment touches the right segment from "above", for mode CA from "below".
 First, "sticking" is assumed, and if this solution does not satisfy the loading condition
 at the contact point, then "sliding" is assumed.
 Admissibility of the solution is checked: 
 * The relative rotation angle at the contact point must be in the appropriate range, otherwise a mode change occurs. 
 * The active length of the segment must not exceed its total length, otherwise a mode change occurs. 
 The return value indicates whether the solution has been found and whether it is admissible. 

 Input variables: 
 * ub_target ... displacements and rotation of the right end wrt to the fixed left end that should be achieved
 * fab ... initial guess of the left-end forces 
 * cmode ... assumed contact mode (must be one of the tip modes)
 * printflag ... flag indicating whether detailed results should be printed

 Output variables: 
 * fab ... left-end forces at the end of the step

 Return value:
 * true or false, indicates whether a converged solution has been found
*/
  ContactModeType findLeftEndForcesLocal_Tip(FloatArray &ub_target, FloatArray &fab, ContactModeType cmode);
  /*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero
and that the process corresponds to ROLLING, i.e., the sum of the active segment lengths
remains constant during the step.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the left end of the beam.
The contact point is located at the inflexion point and its distances from the left and right ends  
(i.e., the active segment lengths) are considered as input as well as output parameters. 
On input, what matters is only their sum. The deformed shape is computed by integrating
over the whole beam because the total length is known. The left-end forces are adjusted
such that the displacements and rotation at the right end are correct. However, if no inflexion
point is detected for the converged solution, the contact conditions cannot be satisfied
and the solution is not admissible. 

Input:
ub_target ... displacements and rotations of the right end with respect to the left end
fab ... initial guess of the left-end forces and moment
deltaPhi ... rotation jump at the contact point (zero for regular modes, +-PI for the inverted modes)
Lac ... initial value of the active length of the left segment
Lbc ... initial value of the active length of the right segment
printflag ... flag indicating whether detailed results should be printed

Output:
fab ... left-end forces and moment at the end of the step
Lac ... active length of the left segment at the end of the step
Lbc ... active length of the right segment at the end of the step
Nc ... normal force at the contact point
Qc ... shear force at the contact point

Return value:
indicates success or failure (failure means that the iterative process does not converge or the deformed
beam does not have an inflexion point) 
*/
bool findLeftEndForcesLocal_Smooth_Rolling(FloatArray &ub_target, FloatArray &fab, double deltaPhi, double &Lac, double &Lbc, double &Nc, double &Qc);
  /*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero
and that the process corresponds to SLIDING, i.e., the sum of the active segment lengths
varies during the step.
Everything is done here in the local coordinate system aligned with the left end of the beam.
The solution is computed iteratively, using the Newton-Raphson method. 
In addition to the left-end forces and moment, the fourth basic unknown is the total beam length
(sum of the active segment lengths). In addition to the conditions dealing with prescribed
displacements and rotation at the right end, the fourth equation is the contact condition
at the contact point. Since sliding is assumed, the combination of the normal force and shear
force must be such that the contact loading function fc vanishes. 

The contact point is located at the inflexion point and its distances from the left and right ends  
(i.e., the active segment lengths) are considered as input as well as output parameters. 
On input, what matters is only their sum. The deformed shape is computed by integrating
over the whole beam because the total length is known. The left-end forces are adjusted
such that the displacements and rotation at the right end are correct. However, if no inflexion
point is detected for the converged solution, the contact conditions cannot be satisfied
and the solution is not admissible. 

Input:
ub_target ... displacements and rotations of the right end with respect to the left end
fab ... initial guess of the left-end forces and moment
deltaPhi ... rotation jump at the contact point (zero for regular modes, +-PI for the inverted modes)
Lac ... initial value of the active length of the left segment
Lbc ... initial value of the active length of the right segment
printflag ... flag indicating whether detailed results should be printed

Output:
fab ... left-end forces and moment at the end of the step
Lac ... active length of the left segment at the end of the step
Lbc ... active length of the right segment at the end of the step
Nc ... normal force at the contact point
Qc ... shear force at the contact point

Return value:
indicates success or failure (failure means that the iterative process does not converge or the deformed
beam does not have an inflexion point) 
*/
bool findLeftEndForcesLocal_Smooth_Sliding(FloatArray &ub_target, FloatArray &fab, double deltaPhi, double &Lac, double &Lbc, ContactModeType cmode);
  /*
Find forces and moment at the left end and the length of the beam that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero, and to the satisfaction
of the conditions describing the sliding unit inside the beam.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the left end of the beam.
The contact point is located at the inflexion point and its distances from the left and right ends  
(i.e., the active segment lengths) are considered as input as well as output parameters. 

Input:
ub_target ... displacements and rotations of the right end with respect to the left end
fab ... initial guess of the left-end forces and moment
cmode ... assumed contact mode (must be one of the smooth ones)
printflag ... flag indicating whether detailed results should be printed

Output:
fab ... left-end forces and moment at the end of the step

Return value:
the same as the assumed mode (input variable cmode) if the solution is admissible, otherwise a guess of the likely contact mode
*/
  ContactModeType findLeftEndForcesLocal_Smooth(FloatArray &ub_target, FloatArray &fab, ContactModeType cmode);
  /*
  This function checks whether the no-contact mode is possible for the given relative displacements.
  If it is, the return value is N_cmode.
  Otherwise the function tries to guess the most likely contact mode but the result is not always
  correct because the behavior is path-dependent.
  The function is used in situation when the mode at the beginning of the step is NOT the no-contact mode.
  Based on the behavior (e.g., lack of convergence or inadmissible solutions for assumed contact modes),
  it is concluded that contact is probably lost, but this assumption needs to be verified, which is
  the purpose of the present function. 
*/
  ContactModeType predictContactMode(FloatArray &ub);
  /*
  This method finds the time at which the tip of a segment hits the x-axis.
  It is assumed that the opposite end of the segment is initially located at (u_prev[0],u_prev[1])
  and rotated by u_prev[2], and that its displacement and rotation linearly increases during the step.
  The 'time' is a dimensionless parameter running from 0 at the beginning to 1 at the end of the step.
  Linear increase of rotation of the opposite end means that the tip travels along a curved trajectory.
  A nonlinear equation describing the condition that the current value of the second coordinate vanishes
  is solved iteratively. 

  Input variables:
  u_prev[3] ... 2 coordinates and rotation at the end of the previous step,
  du[3] ... increments of 2 coordinates and rotation during the step

  Return value:
  dimensionless time at which the tip passes through the x-axis
  (if no converged solution is found, the return value is -1., which is outside the admissible interval) 
 */
  double findTipContactTime(FloatArray &u_prev, FloatArray &du, double segLength);
  /*
  This function checks whether, during the step from ub_prev to ub, contact would occur. 
  If it does, the return value indicates which contact mode would arise, 
  otherwise the result is the N_cmode.
  It is assumed that the state at ub_prev is a no-contact state.
  For contact to occur during the step, the tip of one segment must pass through the other segment.
  The tip that hits the other segment first determines the contact mode.

  Input variables:
  ub[3] ... relative displacements and rotation at the end of the current step
  ub_prev[3] ... relative displacements and rotation at the end of the previous step,

  Return value:
  most likely type of contact mode at the end of the step
*/
    ContactModeType predictContactMode(FloatArray ub, FloatArray ub_prev);
  /*
This method attempts to find (iteratively) the left-end forces and moment that lead to
the given right-end displacements and rotation, taken as relative with respect to the left end.
Everything is done in the local coordinate system attached to the left end.

Input variables:
ub[3] ... relative displacements and rotation at the end of the current step
ub_prev[3] ... relative displacements and rotation at the end of the previous step,
fab[3] ...  initial guess of the left-end forces and moment
printflag ... flag indicating whether detailed results should be printed

Output variables:
fab[3] ...  left-end forces and moment (values of initial guess are rewritten)

Return value:
boolean, indicating whether an admissible solution has been found

The method first assumes that the contact mode remains the
same as in the previous step and attempts to compute the
corresponding solution. If this attempt fails, the assumption
is modified and the whole process is repeated until an admissible
solution is found or until the number of attempts exceeds
a given limit (currently set to 5). If an admissible solution
is found, the key internal variables (the contact mode and the 
active lengths of both segments) are updated, the left-end forces 
at the end of the step are passed as an output variable, and the method
returns Boolean value 'true' to indicate success.

If the assumed mode is no-contact, it is not necessary to solve
a special set of equations, just to check whether the cantilevers
would intersect if they remain straight. If an intersection is 
found, the most likely mode is guessed (based on certain rules),
otherwise the no-contact solution is accepted, which of course
implies that the end forces and element stiffness are set to zero.
The check for intersection of two straight cantilevers is performed 
using the predictContactMode method, which has two different versions,
depending on whether the contact mode at the beginning of the step
is or is not the no-contact mode.
*/
    bool findLeftEndForcesLocal(FloatArray &ub, FloatArray &ub_prev, FloatArray &fab);

    /*
In most cases, the method will be called when the values of fab have already been computed.
However, to make sure that these are really the correct values, they are first recomputed.
If the input contains the correct converged values, the process will now converge after
one iteration. Otherwise, the full iterative process is used and it may not converge.
  */
    void construct_T(FloatMatrix &T, double phia);
  /*
    Auxiliary matrix for transformation of stiffness.
    It is obtained by differentiating T with respect to phia.
  */
  void construct_Tprime(FloatMatrix &T, double phia);
  void construct_l(FloatArray &l, double phia);
  void construct_l_IC(FloatArray &l, double phia, double L);
  void construct_lprime(FloatArray &l, double phia);
  /*
Find forces and moment at the left end that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
Note that the transformation matrix T is affected by angle alpha that specifies the initial beam geometry. 
*/
  bool findLeftEndForces(FloatArray &u, FloatArray &u_prev, FloatArray &fab);
  //void integrateAlongSegmentAndPlot(double fab[3], double Lb, double segmentLength, double u0[2], double T[2][2], FILE* outfile);
  //void plotSegment(double fab[3], double ub[3], bool isLeftSegment, FILE* outfile);
  //void plotResponse(int nstage, int nstep[10], double ustep[3][10], int iplot[10]);

  void computeSegmentDisplacements(FloatMatrix &uMatrix, const FloatArray &fab, double Lb, double segmentLength, const FloatArray &u0, const FloatMatrix &T);
  

  
};
} // end namespace oofem
#endif // nlbem_internalcontact_h
