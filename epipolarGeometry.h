#ifndef EPIPOLARGEOMETRY_H
#define EPIPOLARGEOMETRY_H

#include "basicGeometry.h"
#include <vector>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::JacobiSVD;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;
using std::vector;
using std::cout;
using std::endl;

void DecomposeEssentialMatrix(const cv::Mat& E, const MatrixXd& points0, const MatrixXd& points1, Matrix3d& R, Vector3d& T, MatrixXd& pointsOut, VectorXd& inliers){

    vector<Matrix3d> rotations;
    rotations.push_back(Matrix3d());
    rotations.push_back(Matrix3d());

    std::vector<Vector3d> translations;
    translations.push_back(Vector3d());
    translations.push_back(Vector3d());

    Matrix3d e, U, Vt, W, Z, R1, R2, S;

    W << 0, 1, 0,
        -1, 0, 0,
         0, 0, 1;

    Z << 0, 1, 0,
        -1, 0, 0,
         0, 0, 0;

    cv::cv2eigen(E, e);

    JacobiSVD<Matrix3d> svd(e, ComputeFullU | ComputeFullV);
    U = svd.matrixU();
    Vt = svd.matrixV().transpose();

    R1 = U * W * Vt;
    R2 = U * W.transpose() * Vt;

    if (R1.determinant() < 0)
    {
        R1 *= -1;
    }
    if (R2.determinant() < 0)
    {
        R2 *= -1;
    }

    S = U * Z * U.transpose();

    T(0) = S(2,1);
    T(1) = S(0,2);
    T(2) = S(1,0);

    translations[0] = T;
    translations[1] = -T;

    rotations[0] = R1;
    rotations[1] = R2;

    MatrixXd projMat0(3,4);
    vector<MatrixXd> projMat1;

    projMat0 <<     1,  0,  0,  0,
                    0,  1,  0,  0,
                    0,  0,  1,  0;

    projMat1.push_back(MatrixXd(3,4));
    projMat1.push_back(MatrixXd(3,4));
    projMat1.push_back(MatrixXd(3,4));
    projMat1.push_back(MatrixXd(3,4));

    projMat1[0] << rotations[0], translations[0];
    projMat1[1] << rotations[0], translations[1];
    projMat1[2] << rotations[1], translations[0];
    projMat1[3] << rotations[1], translations[1];

    int bestProjMatSupport = 0, index = -1;

    for(int i=0; i<4; i++){

        MatrixXd points;
        basicGeometry::TriangulatePoints(projMat0, projMat1[i], points0, points1, points);

        MatrixXd proj0 = ( projMat0 * (points.transpose()) ).transpose();
        MatrixXd proj1 = ( projMat1[i] * (points.transpose()) ).transpose();

        VectorXd test0 = proj0.col(2).array() * points.col(3).array();
        VectorXd test1 = proj1.col(2).array() * points.col(3).array();

        VectorXd test(test0.size());

        int numberOfInliers = 0;

        for(int j=0; j<test0.size(); j++){

            test0(j) = (test0(j) > 0)? 1: 0;
            test1(j) = (test1(j) > 0)? 1: 0;

            test(j) = (test0(j)==1 && test1(j)==1)? 1: 0;

            numberOfInliers += test(j);
        }

        if(numberOfInliers > bestProjMatSupport){

            bestProjMatSupport = numberOfInliers;
            index = i;
            inliers = test;

            pointsOut.resize(numberOfInliers, points.cols());

            for(int j=0, k=0; j<inliers.size(), k<numberOfInliers; j++){

                if(inliers(j) == 1){

                    pointsOut.row(k) = points.row(j);
                    k++;
                }
            }
        }
    }

    switch(index){

        case 0:

            R = rotations[0];
            T = translations[0];
            basicGeometry::Homogenize(pointsOut);
            break;

        case 1:

            R = rotations[0];
            T = translations[1];
            basicGeometry::Homogenize(pointsOut);
            break;

        case 2:

            R = rotations[1];
            T = translations[0];
            basicGeometry::Homogenize(pointsOut);
            break;

        case 3:

            R = rotations[1];
            T = translations[1];
            basicGeometry::Homogenize(pointsOut);
            break;

        default:

            R = MatrixXd(3,3);
            R <<    1,  0,  0,
                    0,  1,  0,
                    0,  0,  1;
            T = VectorXd(3);
            T.setZero();
            pointsOut = MatrixXd();
            inliers = VectorXd(points0.rows());
            inliers.setZero();
            break;
    }

    T.normalize();
}

#endif // EPIPOLARGEOMETRY_H
