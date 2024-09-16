#ifndef EXTERNALCONATCTFORCEBODY_H
#define EXTERNALCONATCTFORCEBODY_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"


class externalContactForceBody
{
public:
	externalContactForceBody(elasticPlate &m_plate, timeStepper &m_stepper,
		double m_stiffness, double m_dBar);
	~externalContactForceBody();

	void computeFc();
	void computeJc();

    void setFirstJacobian();

    Vector3d totalContactForce;
    Vector3d totalContactTorque;

    Vector3d bodyPosition;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double d;

    int ind;

    Vector3d dDdEdge;
    Matrix3d Id3;

    double dEnergydD;
    double d2EnergydD2;

    Matrix3d d2DdEdge2;
    Matrix3d d2EdEdge2;

    double stiffness;
    double dBar;

    VectorXi indexArray;
    MatrixXd jacobian;

    Vector3i indexLocal;

    Vector3d x1, x2, x3;
    Vector3d xp;
    Vector3d xProject;

    Vector3d e1, e2, e3, e4;

    Vector3d surfN;

    double dist;

    VectorXd gradientVec;
    MatrixXd hessionMat;

    Matrix3d localJacobian;

    double d1, d2, d3;

    double ifPositive;

    Vector3d pa, pb, pc;

    double minDis;

    double dmin;

    double dLocal;

    int minTriIndex;

    Vector3d xEnd;
    Vector3d xCurrent;

    VectorXd ListVecTP(double a1, double a2, double a3, double a4, double a5, double a6, 
        double a7, double a8, double a9, double a10, double a11, double a12);

    MatrixXd ListMatTP(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4, VectorXd a5, VectorXd a6, 
        VectorXd a7, VectorXd a8, VectorXd a9, VectorXd a10, VectorXd a11, VectorXd a12);

    VectorXd getGradientTP(double xp, double yp, double zp, 
        double x1, double y1, double z1,
        double x2, double y2, double z2,
        double x3, double y3, double z3);

    MatrixXd getHessionTP(double xp, double yp, double zp, 
        double x1, double y1, double z1,
        double x2, double y2, double z2,
        double x3, double y3, double z3);

    VectorXd ListVecEP(double a1, double a2, double a3, double a4, double a5, double a6, 
        double a7, double a8, double a9);

    MatrixXd ListMatEP(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4, VectorXd a5, VectorXd a6, 
        VectorXd a7, VectorXd a8, VectorXd a9);

    VectorXd getGradientEP(double xp, double yp, double zp, 
        double x1, double y1, double z1,
        double x2, double y2, double z2);

    MatrixXd getHessionEP(double xp, double yp, double zp, 
        double x1, double y1, double z1,
        double x2, double y2, double z2);

    VectorXd ListVecEE(double a1, double a2, double a3, double a4, double a5, double a6, 
        double a7, double a8, double a9, double a10, double a11, double a12);

    MatrixXd ListMatEE(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4, VectorXd a5, VectorXd a6, 
        VectorXd a7, VectorXd a8, VectorXd a9, VectorXd a10, VectorXd a11, VectorXd a12);

    double distanceEE(double x11, double y11, double z11,
                           double x12, double y12, double z12, 
                           double x21, double y21, double z21, 
                           double x22, double y22, double z22);

    VectorXd getGradientEE(double x11, double y11, double z11,
                           double x12, double y12, double z12, 
                           double x21, double y21, double z21, 
                           double x22, double y22, double z22);

    MatrixXd getHessionEE(double x11, double y11, double z11,
                          double x12, double y12, double z12, 
                          double x21, double y21, double z21, 
                          double x22, double y22, double z22);

    double ClosestPtSegmentSegment( const Vector3d& p1, const Vector3d& q1, const Vector3d& p2, const Vector3d& q2, double& s, double& t, Vector3d& c1, Vector3d& c2 );
    
    double clamp(double scalarA, double minA, double maxA);

};

#endif
