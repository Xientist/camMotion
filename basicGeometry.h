#ifndef BASICGEOMETRY_H
#define BASICGEOMETRY_H

#include <eigen3/Eigen/Dense>
#include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Vector4d;

namespace basicGeometry {

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

Vector3d EquatorialPointFromQ(Vector4d q){

    double q0, q1, q2, q3, den;

    q0 = q(0);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);

    den = 1 / (1+q3);

    return Vector3d(q0*den, q1*den, q2*den);
}

}

#endif // BASICGEOMETRY_H
