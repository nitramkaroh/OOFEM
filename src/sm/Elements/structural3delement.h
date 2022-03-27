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

#ifndef structural3delement_h
#define structural3delement_h


#include "Elements/nlstructuralelement.h"
#include "fbarelementinterface.h"
#include "Loads/pressurefollowerloadinterface.h"


#define _IFT_Structural3DElement_materialCoordinateSystem "matcs" ///< [optional] Support for material directions based on element orientation.

namespace oofem {
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Base class 3D elements.
 *
 * @author Jim Brouzoulis
 * @author Mikael Öhman
 * @author Martin Horák
 */
 class Structural3DElement : public NLStructuralElement, FbarElementExtensionInterface, PressureFollowerLoadElementInterface
{
protected:
  /** 
   * FEICellGeometry wrapper 
   */
    FEICellGeometry* cellGeometryWrapper;


    bool matRotation;

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Structural3DElement(int n, Domain * d);
    /// Destructor.
    virtual ~Structural3DElement() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialMode giveMaterialMode();

    FEICellGeometry* giveCellGeometryWrapper(TimeStep *tStep = NULL, double alpha = 1);

    virtual int computeNumberOfDofs();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int);
    
    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane);

    void giveMaterialOrientationAt(FloatArray &x, FloatArray &y, FloatArray &z, const FloatArray &lcoords);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveElementParametricCentroid(FloatArray &answer) { answer = {0.0 , 0.0, 0.0};}
    virtual void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType modeType);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, int lowerIndx = 1, int upperIndx = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, double alpha = 0);
    virtual void computeGaussPoints();

     // Edge support
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext)
    { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }     


        
    //virtual IntegrationRule *GetSurfaceIntegrationRule(int); // old
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void giveSurfaceDofMapping(IntArray &answer, int) const;
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);


    // support for pressure follower load interface
    virtual double surfaceEvalVolumeAround(GaussPoint *gp, int iSurf){return this->computeSurfaceVolumeAround(gp, iSurf);}
    virtual void surfaceEvalNmatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp){this->computeSurfaceNMatrixAt(answer, iSurf, gp);}
    virtual void surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdksi,FloatArray &dxdeta, int iSurf, GaussPoint *gp, TimeStep *tStep);
    virtual IntegrationRule* surfaceGiveIntegrationRule(int order, int iSurf) { return this->giveInterpolation()->giveBoundarySurfaceIntegrationRule(order, iSurf);}

    
    
private:
    double dnx(int i, int arg2);
};


} // end namespace oofem
#endif // structural3delement_h
