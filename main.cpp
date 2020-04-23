//#include <QCoreApplication>
#include <QElapsedTimer>
#include <QDebug>

#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>

// #include <Eigen/Dense>        /*(already included in "BasicGeometry.h")*/
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

// #include <opencv2/opencv.hpp>        /*(already included in "BasicGeometry.h")*/
// #include <opencv2/core/eigen.hpp>    /*(already included in "BasicGeometry.h")*/

// #include "basicGeometry.h"           /*(already included in "epipolarGeometry.h")*/
#include "epipolarGeometry.h"
#include "diff.h"

#include "getCorrespondences/GetCorrespondences.cpp"

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

MatrixXd GemanMcClure(VectorXd x){

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

void OptimizeRotationAndTranslation(const VectorXd& vec, const MatrixXd& points0, const MatrixXd& points1, MatrixXd& JREVect, VectorXd& RE, int factor = 1){

    int numEssentialInliers = points0.rows();
    VectorXd p = vec.head(3);
    VectorXd t = vec.tail(vec.size()-3);

    MatrixXd JEVect;
    MatrixXd E = diff::GradientEssentialMatrixwrtVecTrans(p, t, factor, JEVect);

    MatrixXd JREE(numEssentialInliers, 9);
    RE.resize(numEssentialInliers);

    for(int i = 0; i < numEssentialInliers; i++){
        JREE.row(i) = diff::GradientSamsponErrorwrtEssentialMatrix(E, points0.row(i), points1.row(i), RE(i));
    }

    JREVect = JREE*JEVect;
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

        while((file = fopen(filename.c_str(), "r"))){

            i++;
            fclose(file);
            filename = name + "(" + std::to_string(i) +").txt";
        }
    }

    std::ofstream fichier(filename);

    return fichier;
}

// Experimental

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
};
typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

int m_inputs, m_values;

Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

int inputs() const { return m_inputs; }
int values() const { return m_values; }

};

struct my_functor : Functor<double>
{
    my_functor(MatrixXd points0, MatrixXd points1, int factor): Functor<double>(5, points0.rows()), m_points0(points0), m_points1(points1), m_factor(factor) {}
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        MatrixXd JREVect;

        OptimizeRotationAndTranslation(x, m_points0, m_points1, JREVect, fvec, m_factor);

        return 0;
    }

private:
    MatrixXd m_points0, m_points1;
    int m_factor;
};

// ----

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

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

    std::vector<MatrixXd> traj, trajOpt;
    MatrixXd m(4,4);
    m << 1.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;
    traj.push_back(m);
    trajOpt.push_back(m);

    std::cout << "processing..." << std::endl;

    for(int i=0; i<numberOfTests; i++){

//        double time, subtime;
//        QElapsedTimer timer, subtimer;
//        timer.start();

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

        // Estimated Motion/*(already included in "BasicGeometry.h")*/
//        subtimer.start();
        MatrixXd matches = getCorrespondences(imageFile1, imageFile2);
//        subtime = subtimer.elapsed();

        MatrixXd matches1 = matches.block(0, 0, matches.rows(), 2);
        MatrixXd matches2 = matches.block(0, 2, matches.rows(), 2);

        MatrixXd dist = matches1.block(0, 0, matches1.rows(), 1).array().square() + matches1.block(0, 1, matches1.rows(), 1).array().square();

        cv::Mat m1;
        cv::eigen2cv(matches1, m1);
        cv::Mat m2;
        cv::eigen2cv(matches2, m2);
        cv::Mat Kcv(3,3,CV_64FC1);
        cv::eigen2cv(K, Kcv);
        cv::Mat mask;
        cv::Mat E(3,3,CV_64FC1);

        E = cv::findEssentialMat(m1, m2, Kcv, cv::RANSAC, 0.999, 1.0, mask);

//        std::ofstream e("E.txt");

//        for(int k=0; k<3; k++){

//            for(int l=0; l<3; l++){

//                e << E.at<double>(k, l) << " ";
//            }
//        }

//        e.close();

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

        decomposedMatrix dm = DecomposeEssentialMatrix(E, points1, points2);

//        std::cout << "R:" << std::endl;
//        std::cout << dm.R << std::endl;

//        std::cout << "t:" << std::endl;
//        std::cout << dm.t << std::endl;

        // Updating the transformation matrix vector
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

        MatrixXd Pose = traj.back() * transfo;

        traj.push_back(Pose);

        // Optimization
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

        int factor = 1;
        Vector2d v = basicGeometry::EquatorialPointFromT(dm.t, factor);
        VectorXd vec(5);
        vec << u, v;

//        auto fun = [=](VectorXd x) -> VectorXd{
//            MatrixXd JREVect;
//            return OptimizeRotationAndTranslation(x, points1, points2, JREVect, factor);
//        };
//        auto jac = [=](VectorXd x) -> MatrixXd{
//            MatrixXd JREVect;
//            OptimizeRotationAndTranslation(x, points1, points2, JREVect, factor);
//            return JREVect;
//        };
        auto loss = [=](VectorXd x){
            return GemanMcClure(x);
        };

        // Experimental

        my_functor functor(points1, points2, factor);

        Eigen::NumericalDiff<my_functor> numDiff(functor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor> > lm(numDiff);

        lm.parameters.maxfev = 50;
        lm.parameters.ftol = 1.0e-8;
        lm.parameters.gtol = 1.0e-8;
        lm.parameters.xtol = 1.0e-8;

        lm.minimize(vec);

        MatrixXd J;
        MatrixXd ROpt = basicGeometry::Quaternion2Matrix(diff::QFromStereographic(vec, J));
        VectorXd vEnd(2);
        vEnd << vec(3), vec(4);
        VectorXd TOpt = diff::TFromStereographic(vEnd, J, factor);

        VectorXd tOptPoses;
        MatrixXd ROptt = ROpt.transpose();

        tOptPoses = -ROptt * TOpt;
        tOptPoses = tOptPoses * scale;
        tOptPoses.conservativeResize(4);
        tOptPoses(3) = 0;

        MatrixXd transfoOpt = ROptt;
        transfoOpt.conservativeResize(4,4);
        transfoOpt.col(3) = tOptPoses;
        transfoOpt.row(3) << 0, 0, 0, 1;

        MatrixXd PoseOpt = trajOpt.back() * transfoOpt;

        trajOpt.push_back(PoseOpt);

        // ----

        MatrixXd temp;
        VectorXd vec2;
        OptimizeRotationAndTranslation(vec, points1, points2, temp, vec2, factor);
        vec2 = vec2.array() * vec2.array();
        errorStart(i) = 0.5*(loss(vec2)(0));

//        time = timer.elapsed();

//        qDebug() << "Time elapsed for subcode: " << subtime;
//        qDebug() << "Time elapsed for the whole loop: " << time;
//        qDebug() << "Percentage of time taken by subcode: " << subtime / time * 100 << "%";
    }

    std::cout << "Algorithm OK." << std::endl;

    std::ofstream calculated_poses = createFile("calculated_poses");
    std::ofstream optimized_poses = createFile("optimized_poses");

    for(int t=0; t<traj.size(); t++){

        for(int l=0; l<4; l++){

            for(int m=0; m<4; m++){

                calculated_poses << traj[t](l, m) << " ";
                optimized_poses << trajOpt[t](l,m) << " ";
            }
        }

        calculated_poses << std::endl;
        optimized_poses << std::endl;
    }

    calculated_poses.close();
    optimized_poses.close();

    std::cout << "Pose files written." << std::endl;

    errorT = errorT.cwiseAbs();
    errorR = errorR.cwiseAbs();

    //return a.exec();
}
