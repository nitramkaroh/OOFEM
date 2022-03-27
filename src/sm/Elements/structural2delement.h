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

#ifndef structural2delement_h
#define structural2delement_h

#include "Elements/nlstructuralelement.h"
#include "feinterpol2d.h"
#include "fbarelementinterface.h"
#include "Loads/pressurefollowerloadinterface.h"

#define _IFT_Structural2DElement_materialCoordinateSystem "matcs" ///< [optional] Support for material directions based on element orientation.

namespace oofem {

/**
 * Base class for planar 2D elements.
 *
 * @author Jim Brouzoulis
 * @author Martin Horák
 */
  class Structural2DElement : public NLStructuralElement, FbarElementExtensionInterface, PressureFollowerLoadElementInterface
{

protected:
    /** 
     * To facilitate the transformation of 2d elements into 3d, the complexity of transformation from 3d to
     *  local 2d system can be efficiently hidden in custom FEICellGeometry wrapper, that performs the transformation
     * into local system. This way, the existing 2d intrpolation classes can be used.
     * The element maintain its FEICellGeometry, which is accesible through the giveCellGeometryWrapper. 
     * Generalization to 3d then would require only substitution of the geometry warpper and definition of 
     * element transformation matrix.
     */
    FEICellGeometry* cellGeometryWrapper;

    bool matRotation;

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Structural2DElement(int n, Domain * d);
    /// Destructor.
    virtual ~Structural2DElement();
    virtual void postInitialize();
    virtual int giveNumberOfNodes() const;
    /**
     * Returns the Cell Geometry Wrapper. Default inplementation creates FEIElementGeometryWrapper.
     */
    virtual FEICellGeometry* giveCellGeometryWrapper(TimeStep *tStep = NULL, double alpha = 1);
    
    virtual int computeNumberOfDofs();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    
    virtual IRResultType initializeFrom(InputRecord *ir);
    
    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane);
    
    virtual void giveElementParametricCentroid(FloatArray &answer) { answer = {0.0 , 0.0};}


    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType modeType);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    
protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, int lowerIndx = 1, int upperIndx = ALL_STRAINS) = 0 ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, double alpha = 0) = 0;
    virtual void computeGaussPoints();

    void giveMaterialOrientationAt( FloatArray &x, FloatArray &y, const FloatArray &lcoords);
    
    // Edge support
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    virtual void computeEdgeNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);

    Interface * giveInterface(InterfaceType interface);
    // support for pressure follower load interface
    virtual double surfaceEvalVolumeAround(GaussPoint *gp, int iSurf){return this->computeEdgeVolumeAround(gp, iSurf);}
    virtual void surfaceEvalNmatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp){this->computeEdgeNMatrixAt(answer, iSurf, gp);}
    virtual void surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdksi,FloatArray &dxdeta, int iSurf, GaussPoint *gp, TimeStep *tStep);
    void surfaceEvalDeformedNormalAt_fromU(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdksi, int iSurf, GaussPoint *gp, TimeStep *tStep, const FloatArray &vU);
    virtual void surfaceEvalNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdksi, int iSurf, GaussPoint *gp, TimeStep *tStep);

    virtual IntegrationRule* surfaceGiveIntegrationRule(int order, int iSurf) { return this->giveInterpolation()->giveBoundaryEdgeIntegrationRule(order, iSurf);}
    void  surfaceEvalNumericalStiffMatrixAt(FloatMatrix &answer,FloatMatrix &dNdx, FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep);

};



class PlaneStressElement : public Structural2DElement
{
public:
    PlaneStressElement(int n, Domain * d);
    virtual ~PlaneStressElement() { }
    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, double alpha = 0);
};


class PlaneStrainElement : public Structural2DElement
{
public:
    PlaneStrainElement(int n, Domain * d);
    virtual ~PlaneStrainElement() { }
    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, double alpha = 0);
};


class AxisymElement : public Structural2DElement
{
public:
    AxisymElement(int n, Domain * d);
    virtual ~AxisymElement() { }
    virtual MaterialMode giveMaterialMode() { return _3dMat; }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicLength(const FloatArray &crackToNormalPlane);
    virtual double computeVolumeAround(GaussPoint *gp);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep = NULL, double alpha = 0);
    virtual void computeGaussPoints();
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
};


} // end namespace oofem
#endif // structural2delement_h
