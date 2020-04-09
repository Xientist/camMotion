#ifndef EPIPOLARGEOMETRY_H
#define EPIPOLARGEOMETRY_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include "basicGeometry.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::BDCSVD;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;
using std::vector;

struct decomposedMatrix{

    MatrixXd R;
    VectorXd t;
    VectorXd inliers;
    MatrixXd points;
};

decomposedMatrix DecomposeEssentialMatrix(const MatrixXd& E, const MatrixXd& points0, const MatrixXd& points1){

    vector<MatrixXd> rotations;
    rotations.push_back(MatrixXd(3,3));
    rotations.push_back(MatrixXd(3,3));

    std::vector<VectorXd> translations;
    translations.push_back(VectorXd(3));
    translations.push_back(VectorXd(3));

    MatrixXd W(3,3),Z(3,3);

    W <<    0,  1,  0,
           -1,  0,  0,
            0,  0,  1;

    Z <<    0,  1,  0,
           -1,  0,  0,
            0,  0,  0;

    MatrixXd U, S, V;
    BDCSVD<MatrixXd> svd(E, ComputeFullU | ComputeFullV);
    U = svd.matrixU();
    V = svd.matrixV();

    rotations[0] = U*W*V;
    rotations[1] = U*(W.transpose())*V;

    if(rotations[0].determinant() < 0){

        rotations[0] *= -rotations[0];
    }

    if(rotations[1].determinant() < 0){

        rotations[1] *= -rotations[1];
    }

    S = U*Z*(U.transpose());

    VectorXd t(3);
    t << S(2,1), S(0,2), S(1,0);

    translations[0] = t;
    translations[1] = -t;

    MatrixXd projMat0(3,4);
    vector<MatrixXd> projMat1;

    projMat0 <<     1,  0,  0,  0,
                    0,  1,  0,  0,
                    0,  0,  1,  0;MatrixXd(1,1);

    projMat1.push_back(MatrixXd(3,4));
    projMat1.push_back(MatrixXd(3,4));
    projMat1.push_back(MatrixXd(3,4));
    projMat1.push_back(MatrixXd(3,4));

    projMat1[0] << rotations[0], translations[0];
    projMat1[1] << rotations[0], translations[1];
    projMat1[2] << rotations[1], translations[0];
    projMat1[3] << rotations[1], translations[1];

    int bestProjMatSupport = 0, index = -1;
    MatrixXd pointsOut;
    VectorXd inliers;

    for(int i=0; i<4; i++){

        MatrixXd points = basicGeometry::TriangulatePoints(projMat0, projMat1[i], points0, points1);

        MatrixXd proj0 = ( projMat0 * (points.transpose()) ).transpose();
        MatrixXd proj1 = ( projMat1[i] * (points.transpose()) ).transpose();

        VectorXd test0 = proj0.col(2).array() * points.col(3).array();
        VectorXd test1 = proj1.col(2).array() * points.col(3).array();

        VectorXd test(test0.size());

        for(int j=0; j<test0.size(); j++){

            test0(j) = (test0(j) > 0)? 1: 0;
            test1(j) = (test1(j) > 0)? 1: 0;

            test(j) = (test0(j)==1 && test1(j)==1)? 1: 0;
        }

        int numberOfInliers = test.sum();

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

    MatrixXd R, points;

    switch(index){

        case 0:

            R = rotations[0];
            t = translations[0];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        case 1:

            R = rotations[0];
            t = translations[1];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        case 2:

            R = rotations[1];
            t = translations[0];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        case 3:

            R = rotations[1];
            t = translations[1];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        default:

            R = MatrixXd(3,3);
            R <<    1,  0,  0,
                    0,  1,  0,
                    0,  0,  1;
            t = VectorXd(3);
            t.setZero();
            points = MatrixXd();
            inliers = VectorXd(points0.rows());
            inliers.setZero();
            break;
    }

    t.normalize();

    decomposedMatrix dm;
    dm.R = R;
    dm.t = t;
    dm.points = points;
    dm.inliers = inliers;

    return dm;
}

#endif // EPIPOLARGEOMETRY_H
