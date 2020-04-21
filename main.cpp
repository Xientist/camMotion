#include <QCoreApplication>
//#include <opencv2/opencv.hpp>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "getCorrespondences/GetCorrespondences.cpp"
#include <opencv2/core/mat.hpp>
#include <opencv2/core/eigen.hpp>
#include "epipolarGeometry.h"
// #include "basicGeometry.h" (already included in "epipolarGeometry.h")

using std::vector;
using std::string;
using std::ifstream;
using std::istringstream;
//using namespace cv;
using namespace Eigen;

// return data from file in a 2d vector format
// (rows = lines in file, cols = (numeric) elements in line)
vector< vector<double> > loadtxt(string file){

    vector< vector<double> > data;
    ifstream myFile(file);

    string line;

    while(getline(myFile, line)){

        data.push_back(vector<double>());
        istringstream ss(line);

        double value;

        int i=0;
        string temp;

        // fetching all the data on the line
        while(ss >> value){

            data.back().push_back(value);
        }

        // in case a non-numeric value was found in front of the useful data
        if(i==0){

            ss = istringstream(line);

            // fetching the non-numeric data
            ss >> temp;

            // we suppose only numeric data remains on the line
            while(ss >> value){

                data.back().push_back(value);
            }
        }
    }

    return data;
}

MatrixXd GemanMcClure(vector<double> x){

    ArrayXd x1 = Map< ArrayXd >(x.data(), x.size());
    x1 = 1+x1;

    ArrayXd x12(x1.pow(2));

    ArrayXd x13(x1 * x12);

    x1 = 0.5*(x1-1) / x1;
    x12 = 0.5 / x12;
    x13 = -1 / x13;

    ArrayXd temp(x.size()*3);
    temp << x1, x12, x13;

    // Function, First Derivative, Second Derivative
    double *data = temp.data();
    Map< Matrix<double, 3, Dynamic, RowMajor> > output(data, 3, x.size());

    // std::cout << output << std::endl;

    return output;
}

Vector3d rad2deg(Vector3d v){

    double rollDeg, yawDeg, pitchDeg;

    rollDeg = (v(0) * 180) / M_PI;
    yawDeg = (v(1) * 180) / M_PI;
    pitchDeg = (v(2) * 180) / M_PI;

    return Vector3d(rollDeg, yawDeg, pitchDeg);
}

std::ofstream createFile(string name){

    string filename = name + ".txt";

    if(FILE *file = fopen(filename.c_str(), "r")){

        fclose(file);

        int i = 1;
        filename = name + "(" + std::to_string(i) +").txt";

        while(file = fopen(filename.c_str(), "r")){

            i++;
            fclose(file);
            filename = name + "(" + std::to_string(i) +").txt";
        }
    }

    std::ofstream fichier(filename);

    return fichier;
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // Paths to the sequences and poses folder
    int indexSequence = 7;
    string s_index = (indexSequence>9)
            ? std::to_string(indexSequence)
            : "0"+std::to_string(indexSequence);

    string imageFolder = "dataset/sequences/" + s_index + "/image_0/";
    string groundTruthFile = "poses/" + s_index + ".txt";
    string calibFile = "dataset/sequences/" + s_index + "/calib.txt";

    // Loading ground truths from file
    vector< vector<double> > vecPoses;
    vecPoses = loadtxt(groundTruthFile);

    vector<MatrixXd> gt;

    for(int i=0; i<vecPoses.size(); i++){

        Matrix<double, 3, 4> mat = Map< Matrix<double, 3, 4, RowMajor> >(vecPoses[i].data());
        gt.push_back(mat);

        // std::cout << mat << std::endl << std::endl;
    }

    // Initialization of useful matrices
    vector< vector<double> > vec2d;
    vec2d = loadtxt(calibFile);

    double *t2d = vec2d[0].data();

    MatrixXd temp = Map< Matrix<double, 3, 4, RowMajor> >(t2d);
    MatrixXd K = temp.block(0, 0, 3, 3);
    MatrixXd Kinv = K.inverse();

    // std::cout << Kinv << std::endl;

    // IMPORTANT:
    // Streams in C need to have an endline (\n) at the end of a char array to display char data
    // If not, it waits for more data to display until it finds such an endline to display it all at once.
    // Since matrices from eigen are displayed without one after the last row, this row is not displayed
    // Putting an endline manually after displaying a matrix fixes it and displays that last row

    const int numberOfTests = 500;

    // initialization of errors matrices
    MatrixXd errorT     =   Matrix<double, numberOfTests, 3, RowMajor>();  errorT.setZero();
    MatrixXd errorTOpt  =   Matrix<double, numberOfTests, 3, RowMajor>();  errorTOpt.setZero();
    MatrixXd errorR     =   Matrix<double, numberOfTests, 3, RowMajor>();  errorR.setZero();
    MatrixXd errorROpt  =   Matrix<double, numberOfTests, 3, RowMajor>();  errorROpt.setZero();
    MatrixXd errorStart =   Matrix<double, numberOfTests, 1>();  errorStart.setZero();
    MatrixXd errorStop  =   Matrix<double, numberOfTests, 1>();  errorStop.setZero();

    // usage of fastVisualOdometry (python version ln. 91-93)
    // todo

    std::ofstream calculated_poses = createFile("calculated_poses");

    std::vector<MatrixXd> transformations;
    MatrixXd m(4,4);
    m << 1.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;
    transformations.push_back(m);

    for(int i=0; i<numberOfTests; i++){

        int indexImage1 = i;
        int indexImage2 = i+1;

        // paths to the 2 images we are working with in the loop
        string sindex1 = std::to_string(indexImage1);
        string sindex2 = std::to_string(indexImage2);

        const int number_of_zeros = 6;

        sindex1 = std::string(number_of_zeros - sindex1.length(), '0') + sindex1;
        sindex2 = std::string(number_of_zeros - sindex2.length(), '0') + sindex2;

        string imageFile1 = imageFolder + sindex1 + ".png";
        string imageFile2 = imageFolder + sindex2 + ".png";

        // Motion from ground truth
        MatrixXd gt1(4,4);
        gt1 << 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 1.0;
        gt1.row(0) = gt[indexImage1].row(0);
        gt1.row(1) = gt[indexImage1].row(1);
        gt1.row(2) = gt[indexImage1].row(2);

        MatrixXd gt2(4,4);
        gt2 << 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 1.0;
        gt2.row(0) = gt[indexImage2].row(0);
        gt2.row(1) = gt[indexImage2].row(1);
        gt2.row(2) = gt[indexImage2].row(2);

        MatrixXd motionGt = gt1.inverse() * gt2;
        MatrixXd RGt = motionGt.block(0, 0, 3, 3).transpose();
        MatrixXd tGt = -RGt * motionGt.block(0, 3, 3, 1);
        tGt = tGt / tGt.norm();

        double scale = motionGt.block(0, 3, 3, 1).norm();

        // Estimated Motion
        MatrixXd matches = getCorrespondences(imageFile1, imageFile2);
        MatrixXd matches1 = matches.block(0, 0, matches.rows(), 2);
        MatrixXd matches2 = matches.block(0, 2, matches.rows(), 2);

        MatrixXd dist = matches1.block(0, 0, matches1.rows(), 1).array().square() + matches1.block(0, 1, matches1.rows(), 1).array().square();

        cv::Mat m1;
        cv::eigen2cv(matches1, m1);
        cv::Mat m2;
        cv::eigen2cv(matches2, m2);
        cv::Mat Kcv;
        cv::eigen2cv(K, Kcv);
        cv::Mat mask;
        cv::Mat E = cv::findEssentialMat(m1, m2, Kcv, cv::RANSAC, 0.999, 1.0, mask);

        MatrixXd inliers;
        cv::cv2eigen(mask, inliers);

        MatrixXd matchesInliers1(0, 2);
        MatrixXd matchesInliers2(0, 2);
        int l = 0;
        for(int k = 0; k < matches.rows(); k++){
            if(inliers(k) == 1){
                matchesInliers1.conservativeResize(matchesInliers1.rows()+1, 2);
                matchesInliers2.conservativeResize(matchesInliers2.rows()+1, 2);
                matchesInliers1.row(l) = matches1.row(k);
                matchesInliers2.row(l) = matches2.row(k);
                l += 1;
            }
        }

        MatrixXd points1 = (Kinv * (basicGeometry::Homogeneous(matchesInliers1).transpose())).transpose();
        MatrixXd points2 = (Kinv * (basicGeometry::Homogeneous(matchesInliers2).transpose())).transpose();

        MatrixXd Ematrix;
        cv::cv2eigen(E, Ematrix);

        decomposedMatrix dm = DecomposeEssentialMatrix(Ematrix, points1, points2);

        //traj

        int r = points1.rows();
        l = 0;
        MatrixXd temp1(0, points1.cols());
        MatrixXd temp2(0, points2.cols());
        for(int k = 0; k < r; k++){
            if(dm.inliers(k) == 1){
                temp1.conservativeResize(temp1.rows()+1, points1.cols());
                temp2.conservativeResize(temp2.rows()+1, points2.cols());
                temp1.row(l) = points1.row(k);
                temp2.row(l) = points2.row(k);
                l += 1;
            }
        }

        errorT.row(i) = dm.t - tGt;
        errorR.row(i) = rad2deg(basicGeometry::RotationMatrix2PitchYawRoll(dm.R)) - rad2deg(basicGeometry::RotationMatrix2PitchYawRoll(RGt));

        Vector4d q = basicGeometry::Matrix2Quaternion(dm.R);
        Vector3d u = basicGeometry::EquatorialPointFromQ(q);

        VectorXd tPoses;
        MatrixXd Rt = dm.R.transpose();

        tPoses = -Rt * dm.t;
        tPoses = tPoses * scale;
        tPoses.conservativeResize(4);
        tPoses(3) = 0;

        MatrixXd transfo = Rt;
        transfo.conservativeResize(4,4);
        transfo.col(3) = tPoses;
        transfo.row(3) << 0, 0, 0, 1;

        MatrixXd Pose = transformations.back() * transfo;
        transformations.push_back(Pose);
    }

    std::cout << "Algorithm OK." << std::endl;

    for(int t=0; t<transformations.size(); t++){

        for(int l=0; l<4; l++){

            for(int m=0; m<4; m++){

                calculated_poses << transformations[t](l, m) << " ";
            }
        }

        calculated_poses << std::endl;
    }

    calculated_poses.close();

    std::cout << "Pose file written." << std::endl;

    errorT = errorT.cwiseAbs();
    errorR = errorR.cwiseAbs();

    return a.exec();
}
