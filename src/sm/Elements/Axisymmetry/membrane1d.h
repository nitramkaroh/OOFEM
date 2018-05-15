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

#ifndef membrane1d_h
#define membrane1d_h

#include "Elements/nlstructuralelement.h"
#include "Loads/pressurefollowerloadinterface.h"
#include "interfacetype.h"

///@name Input fields for Membrane1d
//@{
#define _IFT_Membrane1d_Name "membrane1d"
#define _IFT_Membrane1d_InitialStretch "initstretch"
//@}

namespace oofem {
class FEI2dLineLin;

/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral -
 * elasticity finite element for axisymmetric 3d continuum.
 * Each node has 2 degrees of freedom.
 */
class Membrane1d : public NLStructuralElement, public PressureFollowerLoadElementInterface
{
protected:
    int pitch;
    double theta;
    double length;
    double initialStretch;
    static FEI2dLineLin interp;
public:
    Membrane1d(int n, Domain * d);
    virtual ~Membrane1d();

    virtual FEInterpolation *giveInterpolation() const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Membrane1d_Name; }
    virtual const char *giveClassName() const { return "Membrane1d"; }

    virtual void computeStiffnessMatrix(FloatMatrix &answer,MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

    virtual int computeNumberOfDofs() { return 4; }
    virtual int computeNumberOfGlobalDofs() { return 4; }
    bool giveRotationMatrix(FloatMatrix &answer);

    double computeMidsurfaceVolume(TimeStep *tStep);


    
    IRResultType initializeFrom(InputRecord *ir);
protected:

    void computeNuMatrixAt(const FloatArray &iLocCoord, FloatArray &answer);
    void computeNwMatrixAt(const FloatArray &iLocCoord, FloatArray &answer);
    void computeBuMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep);
    void computeBwMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep);



    void computeBrMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep, bool linearized);
    void computeBtMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep, bool linearized);


    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool linearized = false) ;
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int lowerIndx = 1, int upperIndx = ALL_STRAINS){;}
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    void computeG1G2MatricesAt(FloatMatrix &G1, FloatMatrix &G2, GaussPoint *gp, TimeStep *tStep);

    void postInitialize();
    void computeGaussPoints();
    virtual MaterialMode giveMaterialMode() { return _AxisymMemebrane1d; }
    double computeLength();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);

    // edge load support
    void resolveCoordIndices(int &c1, int &c2);
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void giveEdgeDofMapping(IntArray &answer, int) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *gp);
    double givePitch();
    bool computeGtoLRotationMatrix(FloatMatrix &answer);
    virtual Interface *giveInterface(InterfaceType it);


    // pressure follower load load support

    virtual double surfaceEvalVolumeAround(GaussPoint *gp, int iSurf){ return this->computeVolumeAround(gp); }
    virtual void surfaceEvalNmatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp){ this->computeNmatrixAt(gp->giveNaturalCoordinates(), answer); }
    virtual void surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdksi, FloatArray &dxdeta, int iSurf, GaussPoint *gp, TimeStep *tStep);
    virtual IntegrationRule* surfaceGiveIntegrationRule(int order, int iSurf){ return this->giveDefaultIntegrationRulePtr(); }
    virtual void surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp){;}

    

};
} // end namespace oofem
#endif // membrane1d
