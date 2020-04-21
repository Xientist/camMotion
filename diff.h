#ifndef DIFF_H
#define DIFF_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include "basicGeometry.h"
#include <vector>
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd JacobianRotationwrtQuaternion(VectorXd q){
    MatrixXd J(9, 4);
    double w2 = 2.0*q[0];
    double x2 = 2.0*q[1];
    double y2 = 2.0*q[2];
    double z2 = 2.0*q[3];

    double x4 = 2.0*x2;
    double y4 = 2.0*y2;
    double z4 = 2.0*z2;

    J << 0, 0, -y4, -z4,
        -z2, y2, x2, -w2,
         y2, z2, w2, x2,
         z2, y2, x2, w2,
         0, -x4, 0, -z4,
        -x2, -w2, z2, y2,
        -y2, z2, -w2, x2,
         x2, w2, z2, y2,
         0, -x4, -y4, 0;

    return J;
}

VectorXd GradientSamsponErrorwrtEssentialMatrix(MatrixXd E, VectorXd pointImage1, VectorXd pointImage2, double &error){
    double e1 = E(0, 0);
    double e2 = E(0, 1);
    double e3 = E(0, 2);
    double e4 = E(1, 0);
    double e5 = E(1, 1);
    double e6 = E(1, 2);
    double e7 = E(2, 0);
    double e8 = E(2, 1);
    double e9 = E(2, 2);

    double x = pointImage1[0];
    double y = pointImage1[1];
    double xp = pointImage2[0];
    double yp = pointImage2[1];

    double num1 = e9 + e3*xp + e6*yp;
    double num2 = e7 + e1*xp + e4*yp;
    double num3 = e8 + e2*xp + e5*yp;

    double num = num1 + x*num2 + y*num3;

    double den1 = e3 + e1*x + e2*y;
    double den2 = e6 + e4*x + e5*y;
    double den3 = e7 + e1*xp + e4*yp;
    double den4 = e8 + e2*xp + e5*yp;

    double denT = (den1*den1) + (den2*den2) + (den3*den3) + (den4*den4);

    double numdenT32 = num * pow(denT, -1.5);

    double den = 1.0/sqrt(denT);

    error = num*den;

    VectorXd gradient(9);
    gradient[0] = den*x*xp - numdenT32 * (x*den1+xp*den3);
    gradient[1] = den*y*xp - numdenT32 * (y*den1+xp*den4);
    gradient[2] = den*xp - numdenT32 * den1;
    gradient[3] = den*yp*x - numdenT32 * (yp*den3+x*den2);
    gradient[4] = den*yp*y - numdenT32 * (y*den2+yp*den4);
    gradient[5] = den*yp - numdenT32 * den2;
    gradient[6] = den*x - numdenT32 * den3;
    gradient[7] = den*y - numdenT32 * den4;
    gradient[8] = den;

    return gradient;
 }

MatrixXd GradientEssentialMatrixwrtVecTrans(VectorXd point, VectorXd t, int factor, MatirxXd &JEVect){
    MatrixXd JqVec;
    VectorXd q = QFromStereographic(point, JqVec);
    MatrixXd R = basicGeometry.Quaternion2Matrix(q);
    double r1 = R(0, 0);
    double r2 = R(0, 1);
    double r3 = R(0, 2);
    double r4 = R(1, 0);
    double r5 = R(1, 1);
    double r6 = R(1, 2);
    double r7 = R(2, 0);
    double r8 = R(2, 1);
    double r9 = R(2, 2);

    MatrixXd Jtp;
    t = TFromStereographic(t, Jtp, factor);
    double a = t[0];
    double b = t[1];
    double c = t[2];

    MatrixXd JER = MatrixXd::Zero(9, 9);
    JER(0, 3) = -c;
    JER(0, 6) = b;
    JER(1, 7) = b;
    JER(1, 4) = -c;
    JER(2, 8) = b;
    JER(2, 5) = -c;
    JER(3, 0) = c;
    JER(3, 6) = -a;
    JER(4, 1) = c;
    JER(4, 7) = -a;
    JER(5, 2) = c;
    JER(5, 8) = -a;
    JER(6, 3) = a;
    JER(6, 0) = -b;
    JER(7, 4) = a;
    JER(7, 2) = -b;
    JER(8, 5) = a;
    JER(8, 2) = -b;

    MatrixXd JRq = JacobianRotationwrtQuaternion(q);
    MatrixXd JRVec = JRq*JqVec;

    MatrixXd JEt(9, 3);
    JEt << 0, r7, -r4,
           0, r8, -r5,
           0, r9, -r6,
           -r7, 0, r1,
           -r8, 0, r2,
           -r9, 0, r3,
           r4, -r1, 0,
           r5, -r2, 0,
           r6, -r3, 0;

    MatrixXd JEp = JEt*Jtp;

    MatrixXd JEVec = JER*JRVec;

    JEVect = MatrixXd(JEVec.rows(), JEVec.cols()+JEp.cols());
    JEVect << JEVec, JEp;

    MatrixXd E = BasicGeometry.CrossProductMatrix(t)*R;

    return E;
}

VectorXd QFromStereographic(VectorXd point, MatrixXd &J){
    double x = point[0];
    double y = point[1];
    double z = point[2];

    double alpha2 = (x*x) + (y*y) + (z*z);

    double den = 1.0/(alpha2 + 1.0);

    VectorXd<double> q(4);
    q << 2.0*x*den, 2.0*y*den, 2.0*z*den, (1-alpha2)*den;

    J = MatrixXd<double>(4, 3);

    double alpha21 = alpha2 + 1.0;
    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    double factor = -4/(alpha21*alpha21);

    J << 0.5*((2.0*(x*x))-alpha21), xy, xz, x,
         xy, 0.5*((2.0*(y*y))-alpha21), yz, y,
         xz, yz, 0.5*((2.0*(z*z))-alpha21), z;

    J = J*factor;

    return q;
}

VectorXd TFromStereographic(VectorXd point, MatrixXd &J, int d=1){
    double x = point[0];
    double y = point[1];

    double alpha2 = (x*x) + (y*y);
    double den = (d*1.0)/(1.0 + alpha2);

    VectorXd<double> t(3);
    t << (1.0-alpha2)*den, 2*x*den, 2*y*den;

    J = MatrixXd<double>(3, 2);

    double alpha21 = alpha2 + 1.0;
    xy = x*y;
    double factor = (-4.0*d)/(alpha21*alpha21);

    J << x, 0.5*((2*(x*x))-alpha21), xy,
         y, xy, 0.5*((2.0*(y*y))-alpha21);

    J = J*factor;

    return t;
}

#endif // DIFF_H
