#ifndef BASICGEOMETRY_H
#define BASICGEOMETRY_H

#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;
using cv::eigen2cv;

namespace basicGeometry {

MatrixXd CrossProductMatrix(VectorXd vec){
    MatrixXd crossMat = MatrixXd::Zero(3, 3);
    //1st line
    crossMat(0, 1) = -vec[2];
    crossMat(0, 2) = vec[1];
    //2nd line
    crossMat(1, 0) = vec[2];
    crossMat(1, 2) = -vec[0];
    //3rd line
    crossMat(2, 0) = -vec[1];
    crossMat(2, 1) = vec[0];

    return crossMat;
}

MatrixXd Homogeneous(const MatrixXd& M){

    MatrixXd mat = M;
    mat.conservativeResize(mat.rows(), mat.cols()+1);

    VectorXd vec = VectorXd(mat.rows());
    vec.setOnes();

    mat.col(mat.cols()-1) = vec;

    return mat;
}

Vector3d RotationMatrix2PitchYawRoll(const MatrixXd& R){

    double roll, yaw, pitch;

    roll    = atan2( R(2,1), R(1,1) );
    yaw     = atan2( -R(2,0), sqrt( R(2,1)*R(2,1) + R(2,2)*R(2,2) ) );
    pitch   = atan2( R(2,1), R(2,2) );

    return Vector3d(roll, yaw, pitch);
}

MatrixXd Quaternion2Matrix(VectorXd q){
    MatrixXd M(3, 3);
    double w = q[0];
    double x = q[1];
    double y = q[2];
    double z = q[3];

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;

    double xy = x*y;
    double xz = x*z;
    double xw = x*w;
    double yz = y*z;
    double yw = y*w;
    double zw = z*w;

    M << 1.0-2.0*(y2+z2), 2.0*(xy-zw), 2.0*(xz+yw),
         2.0*(xy+zw), 1.0-2.0*(x2+z2), 2.0*(yz-xw),
         2.0*(xz-yw), 2.0*(yz+xw), 1.0-2.0*(x2+y2);

    return M;
}

Vector4d Matrix2Quaternion(const MatrixXd& M){

    Vector4d q;
    double t;

    if (M(2,2) < 0){

            if (M(0,0) > M(1,1)){

                // std::cout << "Case 1" << std::endl;
                t = 1 + M(0,0) - M(1,1) - M(2,2);
                q = Vector4d(t, M(0,1)+M(1,0), M(2,0)+M(0,2), M(2,1)-M(1,2));
            } else {

                // std::cout << "Case 2" << std::endl;
                t = 1 - M(0,0) + M(1,1) - M(2,2);
                q = Vector4d(M(0,1)+M(1,0), t, M(1,2)+M(2,1), M(0,2)-M(2,0));
            }

    } else {

            if (M(0,0) < -M(1,1)){

                // std::cout << "Case 3" << std::endl;
                t = 1 - M(0,0) - M(1,1) + M(2,2);
                q = Vector4d(M(2,0)+M(0,2), M(1,2)+M(2,1), t, M(1,0)-M(0,1));

            } else {

                // std::cout << "Case 4" << std::endl;
                t = 1 + M(0,0) + M(1,1) + M(2,2);
                q = Vector4d(M(2,1)-M(1,2), M(0,2)-M(2,0), M(1,0)-M(0,1), t);
            }
    }

    q = Vector4d(q(3), q(0), q(1), q(2));
    q *= 0.5 / sqrt(t);

    return q;
}

Vector3d EquatorialPointFromQ(const Vector4d& q){

    double q0, q1, q2, q3, den;

    q0 = q(0);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);

    den = 1 / (1+q3);

    return Vector3d(q0*den, q1*den, q2*den);
}

Vector2d EquatorialPointFromT(Vector3d t, int d=1){
    double t0 = t(0);
    double t1 = t(1);
    double t2 = t(2);

    double factor = 1.0/(1.0 + t0);

    return Vector2d(t1*factor, t2*factor);
}

MatrixXd TriangulatePoints(const MatrixXd& projMat0, const MatrixXd& projMat1, const MatrixXd& points0, const MatrixXd& points1){

    MatrixXd points(points0.rows(), 4);

    MatrixXd A(4,4);
    A.setZero();

    for(int i=0; i<points0.rows(); i++){

        A.row(0) = points0(i, 1) * projMat0.row(2) - projMat0.row(1);
        A.row(1) = - points0(i, 0) * projMat0.row(2) + projMat0.row(0);
        A.row(2) = points1(i, 1) * projMat1.row(2) - projMat1.row(1);
        A.row(3) = - points1(i, 0) * projMat1.row(2) + projMat1.row(0);

        cv::Mat a(4,4,CV_64FC1), vt(4,4,CV_64FC1);
        cv::Mat w,u;
        cv::eigen2cv(A, a);

        cv::SVDecomp(a,w,u,vt);

        MatrixXd V;
        cv::cv2eigen(vt.t(), V);

        points.row(i) = V.col(V.cols()-1);
    }

    return points;
}

MatrixXd Homogenize(const MatrixXd& M){

    VectorXd temp = M.col(M.cols()-1);
    MatrixXd result = M;

    for(int i=0; i< M.rows(); i++){

        result.row(i) = result.row(i) / temp(i);
    }

    return result;
}

}

#endif // BASICGEOMETRY_H
