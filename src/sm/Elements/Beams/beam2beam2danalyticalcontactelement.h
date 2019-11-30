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

#ifndef beam2beam2danalyticalcontactelement_h
#define beam2beam2danalyticalcontactelement_h

#include "../sm/Elements/Beams/beam2d.h"

///@name Input fields for Beam2d
//@{
#define _IFT_Beam2Beam2dAnalyticalContactElement_Name "beam2beamcontact"
#define _IFT_Beam2Beam2dAnalyticalContactElement_b1 "b1"
#define _IFT_Beam2Beam2dAnalyticalContactElement_b2 "b2"

//@}

namespace oofem {

/**
 * description of the class, i.e., element
 */
class Beam2Beam2dAnalyticalContactElement : public Beam2d
{
protected:
  FloatArray b1;
  FloatArray b2;
  FloatArray n1;
  FloatArray n2;
  
  double lengthBeam1;
  double lengthBeam2;
  double totalLength;
  double overlap;
  double E = -1;
  double Iy;


    /** Type characterizing the contact case
     */
    enum ContactType {
        CT_None = 0,
	CT_Beam1_Tip_Beam2 = 1,
	CT_Beam2_Tip_Beam1 = 2,
	CT_Overlap = 3,
    };


  
public:
    Beam2Beam2dAnalyticalContactElement(int n, Domain *aDomain);
    virtual ~Beam2Beam2dAnalyticalContactElement(){;}

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    void  postInitialize() override;
    IRResultType initializeFrom(InputRecord *ir) override;
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    void updateYourself(TimeStep *tStep);
 protected:
    void checkContact(double l1, double l2, double o, double &xc, TimeStep *tStep, ContactType &ct);
    bool checkOverlapingContact(double l1, double l2, double overlap, double &xc, TimeStep *tStep);
    void checkTipContact(double &answer, double l1, double l2, TimeStep *tStep);

    void giveInternalForcesVector_Contact(FloatArray &answer, double l1, double l2, TimeStep *tStep);
    void computeStiffnessMatrix_TipContact(FloatMatrix &answer, double l1, double l2);
    void computeStiffnessMatrix_OverlapingContact(FloatMatrix &answer, double L);
    void giveLocalDofs(double &w1, double &phi1, double &w2, double &phi2, TimeStep *tStep);
    void printOutputAt(FILE *File, TimeStep *tStep);
    double initYoungModulus( TimeStep *tStep);
    

};
} // end namespace oofem
#endif // beam2beam2danalyticalcontactelement_h
