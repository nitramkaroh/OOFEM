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

#ifndef nlbeam_mod_h
#define nlbeam_mod_h


#include "dofmanager.h"
#include "Elements/nlstructuralelement.h"
#include "scalarfunction.h"
#include "vtkxmlexportmodule.h"


///@name Input fields for NlBeam_SM2
//@{
#define _IFT_NlBeamMod_Name "nlbeam_mod"
#define _IFT_NlBeamMod_NIP "nip"
#define _IFT_NlBeamMod_EA "ea"
#define _IFT_NlBeamMod_EI "ei"
#define _IFT_NlBeamMod_GAs "gas"
#define _IFT_NlBeamMod_initialCurvature "initialcurvature"
#define _IFT_NlBeamMod_RADIUS "radius"
#define _IFT_NlBeamMod_DEPTH "depth"
#define _IFT_NlBeamMod_Material "materialtype"
#define _IFT_NlBeamMod_Tolerance "btol"
#define _IFT_NlBeamMod_MaxIteration "bmaxit"
#define _IFT_NlBeamMod_NumberMaxSubsteps "nsubsteps"
#define _IFT_NlBeamMod_Section_Tolerance "stol"
#define _IFT_NlBeamMod_Section_MaxIteration "smaxit"

#define _IFT_NlBeamMod_pressure "pressure"
#define _IFT_NlBeamMod_pressureLTF "pressure_ltf"

#define _IFT_NlBeamMod_s "s"
#define _IFT_NlBeamMod_u0 "u0"
#define _IFT_NlBeamMod_w0 "w0"
#define _IFT_NlBeamMod_phi0 "phi0"
#define _IFT_NlBeamMod_kappa0 "kappa0"
#define _IFT_NlBeamMod_tangentVector "tv"

#define _IFT_NlBeamMod_beta "beta"
#define _IFT_NlBeamMod_curvedbeamLength "curvedbeamlength"
#define _IFT_NlBeamMod_coordinateFlag "cf"


#define _IFT_NlBeamMod_coupling "coupling"

//@}

namespace oofem {


/**
 * This class implements a 2-dimensional large strain beam element
 * The shooting method is used to calculate internal forces and stiffness matrix
 * Add more description 
 */
  class NlBeamMod : public NLStructuralElement, public VTKXMLExportModuleElementInterface
{
protected:
    int NIP = 100;
    double pitch = 10, beamLength = 0;
    FloatArray internalForces, fab_init;
    FloatArray s, ds, u, w, phi, kappa;
    FloatMatrix jacobi;
    double beam_tol = 1.e-6, beam_maxit = 100;
    double section_tol = 1.e-6,section_maxit = 20;
    int nsubsteps_init = 4;

  double EI, EA, GAs;
    double RADIUS, DEPTH;
    double curvedbeamLength;
    double choordLength;
    int materialtype = 1;
    FloatArray vM, vV, vN;
    double cosBeta, sinBeta, cosAlpha, sinAlpha;
    /// Initial stress field curved configuration
    ScalarFunction u_0, w_0, phi_0, kappa_0, sx;
  //    std::vector<ScalarFunction> u_0, w_0, phi_0, kappa_0, sx;
    FloatArray u0, w0, phi0, phi0mid, kappa0;
    double beta = 0;
    FloatArray tangentVector;
    enum CoordinateFlag{CF_s = 0, CF_x = 1};
    CoordinateFlag cf;
    int coupling = 1;
    double pressure = 0;
    int pressure_ltf;
    
public:
    NlBeamMod(int n, Domain *aDomain);
    virtual ~NlBeamMod(){;}

    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    virtual int computeNumberOfDofs() { return 6; }

    virtual void  printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "NlBeamMod"; }
    virtual const char *giveInputRecordName() const { return _IFT_NlBeamMod_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    void giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
    void giveCompositeExportData_curved(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );
    void giveCompositeExportData_straight(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );

     Interface *giveInterface(InterfaceType it) override;
     Element_Geometry_Type giveGeometryType() const override { return EGT_Composite; }

     void updateYourself(TimeStep *tStep) override;
     void initForNewStep() override;


protected:

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
			  TimeStep *tStep = NULL, int lowerIndx = 1, int upperIndx = ALL_STRAINS){;}
    virtual double computeLength();
    double givePitch();
    virtual MaterialMode giveMaterialMode() { return _2dBeam; }

    double computeMomentFromCurvature(double kappa);
    double computeDerMomentFromCurvature(double kappa);
    double computeCurvatureFromMoment(double M);

    double computeChi_mid(double Xab, double Zab); 
  
    void integrateAlongBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi, TimeStep *tStep);
    void integrateAlongStraightBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi, TimeStep *tStep);
    void integrateAlongCurvedBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi);
    bool findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc, TimeStep *tStep);
    void construct_T(FloatMatrix &T, const double phia);
    void construct_Tprime(FloatMatrix &T, const double phia);
    void construct_l(FloatArray &l, double phia);
    void construct_l(FloatArray &l, double phia, double L);
    void construct_l(FloatArray &l, double phia, double L, double cB, double sB);
 
    void construct_lprime(FloatArray &l, const double phia);
  
    void findLeftEndForces(const FloatArray &u, FloatArray &fab, TimeStep *tStep);

    void  printOutputAt_StraightBeam(FILE *file, TimeStep *tStep);
    void  printOutputAt_CurvedBeam(FILE *file, TimeStep *tStep);

    double eval_kappa0(double x);
    double eval_phi0(double x);
    double eval_u0(double x);
    double eval_w0(double x);
    double eval_s(double x);

    
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

    void findInitialShape(FloatArray &ub_target, FloatArray &fab_loc, const FloatArray &ds, const FloatArray &kappa_0, FloatArray &u_0, FloatArray &w_0, FloatArray &phi_0, double NIP_0);
    void integrateAlongBeamAndFindInitialShape(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi, const FloatArray &ds, FloatArray &kappa_0, FloatArray &u_0, FloatArray &w_0, FloatArray &phi_0, double NIP_0);

    
    void computeStiffnessMatrix_num(FloatMatrix &answer,FloatMatrix &answer_p, MatResponseMode rMode, TimeStep *tStep);    
    void giveInternalForcesVector_from_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u);
    void giveInternalForcesVectorPressure_from_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u);


};
} // end namespace oofem
#endif // beam2d_h
