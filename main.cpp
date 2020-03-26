#include <QCoreApplication>
//#include <opencv2/opencv.hpp>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <eigen3/Eigen/Dense>

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

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // Paths to the sequences and poses folder
    int indexSequence = 0;
    string s_index = (indexSequence>9)
            ? std::to_string(indexSequence)
            : "0"+std::to_string(indexSequence);

    string imageFolder = "dataset/sequences/" + s_index + "/image0/";
    string groundTruthFile = "poses/" + s_index + ".txt";
    string calibFile = "dataset/sequences/" + s_index + "/calib.txt";

    // Loading ground truths from file
    vector< vector<double> > vecPoses;
    vecPoses = loadtxt(groundTruthFile);

    vector<MatrixXd> gt;

    for(int i=0; i<vecPoses.size(); i++){

        // might need to remove RowMajor?
        Matrix<double, 3, 4> mat = Map< Matrix<double, 3, 4, RowMajor> >(vecPoses[i].data());
        gt.push_back(mat);

        // std::cout << mat << std::endl << std::endl;
    }

    // Initialization of useful matrices
    vector< vector<double> > vec2d;
    vec2d = loadtxt(calibFile);

    double *t2d = vec2d[0].data();

    // might need to remove RowMajor?
    MatrixXd temp = Map< Matrix<double, 3, 4, RowMajor> >(t2d);
    MatrixXd K = temp.block(0, 0, 3, 3);
    MatrixXd Kinv = K.inverse();

    // std::cout << Kinv << std::endl;

    // IMPORTANT:
    // Streams in C need to have an endline (\n) at the end of a char array to display char data
    // If not, it waits for more data to display until it finds such an endline to display it all at once.
    // Since matrices from eigen are displayed without one after the last row, this row is not displayed
    // Putting an endline manually after displaying a matrix fixes it and displays that last row

    int numberOfTests = 100;

    // initialization of errors matrices
    // todo

    for(int i=0; i<numberOfTests; i++){

        int indexImage1 = i;
        int indexImage2 = i+1;

        // todo
    }

    return a.exec();
}