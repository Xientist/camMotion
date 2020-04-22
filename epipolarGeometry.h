#ifndef EPIPOLARGEOMETRY_H
#define EPIPOLARGEOMETRY_H

#include "basicGeometry.h"
#include <vector>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::VectorXd;
using Eigen::JacobiSVD;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;
using std::vector;
using std::cout;
using std::endl;

struct decomposedMatrix{

    MatrixXd R;
    VectorXd t;
    VectorXd inliers;
    MatrixXd points;
};

decomposedMatrix DecomposeEssentialMatrix(const cv::Mat& E, const MatrixXd& points0, const MatrixXd& points1){

    vector<MatrixXd> rotations;
    rotations.push_back(MatrixXd(3,3));
    rotations.push_back(MatrixXd(3,3));

    std::vector<VectorXd> translations;
    translations.push_back(VectorXd(3));
    translations.push_back(VectorXd(3));

//    JacobiSVD<MatrixXd> svd(E, ComputeFullU | ComputeFullV);

//    U = svd.matrixU();
//    S = svd.singularValues();
//    V = svd.matrixV();

//    cv::Mat E(3,3,CV_64FC1);
    cv::Mat W(3,3,CV_64FC1, double(0)), Z(3,3,CV_64FC1, double(0));
    cv::Mat R1(3,3,CV_64FC1), R2(3,3,CV_64FC1), S(3,3,CV_64FC1);

//    std::string filename = "E.txt";
//    std::ifstream fileStream(filename);

//    int cols = 3;
//    double m;
//    int cnt = 0;//index starts from 0
//    while (fileStream >> m)
//    {
//        int temprow = cnt / cols;
//        int tempcol = cnt % cols;
//        E.at<double>(temprow, tempcol) = m;
//        cnt++;
//    }

//    fileStream.close();

    cv::Mat t(3,1,CV_64FC1);

    W.at<double>(0,1) = 1;
    W.at<double>(1,0) = -1;
    W.at<double>(2,2) = 1;

    Z.at<double>(1,0) = -1;
    Z.at<double>(0,1) = 1;

    cv::Mat w,u,vt;

    cv::SVDecomp(E,w,u,vt);

//    std::cout << w << std::endl;
//    std::cout << u << std::endl;
//    std::cout << vt << std::endl;

    R1 = u * W * vt;
    R2 = u * W.t() * vt;

    if (cv::determinant(R1) < 0)
    {
        R1 *= -1;
    }
    if (cv::determinant(R2) < 0)
    {
        R2 *= -1;
    }

    S = u * Z * u.t();

//    std::cout << "R1" << std::endl;
//    std::cout << R1 << std::endl;

//    std::cout << "R2" << std::endl;
//    std::cout << R2 << std::endl;

    S = u * Z * u.t();
//    std::cout << "S" << std::endl;
//    std::cout << S << std::endl;

    t.at<double>(0,0) = S.at<double>(2,1);
    t.at<double>(1,0) = S.at<double>(0,2);
    t.at<double>(2,0) = S.at<double>(1,0);

//    std::cout << "t" << std::endl;
//    std::cout << t << std::endl;

    VectorXd T;
    cv::cv2eigen(t, T);

    translations[0] = T;
    translations[1] = -T;

    cv::cv2eigen(R1, rotations[0]);
    cv::cv2eigen(R2, rotations[1]);

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
    MatrixXd pointsOut;
    VectorXd inliers;

    for(int i=0; i<4; i++){

        MatrixXd points = basicGeometry::TriangulatePoints(projMat0, projMat1[i], points0, points1);

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

    MatrixXd R, points;

    switch(index){

        case 0:

            R = rotations[0];
            T = translations[0];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        case 1:

            R = rotations[0];
            T = translations[1];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        case 2:

            R = rotations[1];
            T = translations[0];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        case 3:

            R = rotations[1];
            T = translations[1];
            points = basicGeometry::Homogenize(pointsOut);
            break;

        default:

            R = MatrixXd(3,3);
            R <<    1,  0,  0,
                    0,  1,  0,
                    0,  0,  1;
            T = VectorXd(3);
            T.setZero();
            points = MatrixXd();
            inliers = VectorXd(points0.rows());
            inliers.setZero();
            break;
    }

    T.normalize();

    decomposedMatrix dm;
    dm.R = R;
    dm.t = T;
    dm.points = points;
    dm.inliers = inliers;

    return dm;
}

#endif // EPIPOLARGEOMETRY_H
