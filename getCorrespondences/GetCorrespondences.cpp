#include "isprocessor2.h"
#include "opencv2/opencv.hpp"
#include "string.h"
#include "iostream"
#include <Eigen/Dense>

using namespace Eigen;

#define CORNERS_DATA 1024
#define CDIMAGE_MEANWIDTH 7


MatrixXd getCorrespondences(std::string fileImage1, std::string fileImage2)
{

//    std::cout << fileImage1 << std::endl;
//    std::cout << fileImage2 << std::endl;

    cv::Mat image1OpenCV = cv::imread(fileImage1, cv::IMREAD_GRAYSCALE);
    cv::Mat image2OpenCV = cv::imread(fileImage2, cv::IMREAD_GRAYSCALE);

    /*
    cv::Size size1 = image1OpenCV.size();
    size1.width = size1.width / 2.0;
    size1.height = size1.height / 2.0;

    cv::Size size2 = image2OpenCV.size();
    size2.width = size2.width / 2.0;
    size2.height = size2.height / 2.0;

    cv::resize(image1OpenCV, image1OpenCV, size1);
    cv::resize(image2OpenCV, image2OpenCV, size2);
    */

    int width(image1OpenCV.cols), height(image1OpenCV.rows);

//    std::cout << "width x height " << width << "x" << height << std::endl;

    unsigned char *image1 = new unsigned char[width * height];
    unsigned char *image2 = new unsigned char[width * height];

    unsigned char *it_image1 = image1;
    unsigned char *it_image2 = image2;

    cv::MatIterator_<unsigned char> it, end;

    for (it = image1OpenCV.begin<unsigned char>(), end = image1OpenCV.end<unsigned char>(); it != end; it++)
    {
        *it_image1++ = *it;
    }
    for (it = image2OpenCV.begin<unsigned char>(), end = image2OpenCV.end<unsigned char>(); it != end; it++)
    {
        *it_image2++ = *it;
    }

    // Process
    ISProcessor2* isprocessor = new ISProcessor2(width, height, CORNERS_DATA);

    // Find corners image 1
    IntegralImage integralImage1(width, height);
    integralImage1.setImage(image1);
    integralImage1.meanFilter(CDIMAGE_MEANWIDTH);

    CDCornerData* cornerImage1 = new CDCornerData[CORNERS_DATA];
    isprocessor->findCorners(image1, cornerImage1, integralImage1.imageBuffer());
    
    // Find corners image 2
    IntegralImage integralImage2(width, height);
    integralImage2.setImage(image2);
    integralImage2.meanFilter(CDIMAGE_MEANWIDTH);
    
    CDCornerData* cornerImage2 = new CDCornerData[CORNERS_DATA];
    isprocessor->findCorners(image2, cornerImage2, integralImage2.imageBuffer());

    // Corner matching
    isprocessor->setImage1(image1, cornerImage1);
    isprocessor->setCornersImage1(cornerImage1);
    
    isprocessor->setImage2(image2, cornerImage2);
    isprocessor->setCornersImage2(cornerImage2);

    isprocessor->process();

//    std::cout << "FB corner " << isprocessor->numberOfFBCorners() << std::endl;

    ISProcessor2::Correspondence* correspondences = isprocessor->getCorrespondences();

    int numberOfCorrespondences = isprocessor->getNumberOfCorrespondences();

    //std::cout << numberOfCorrespondences << std::endl;
    MatrixXd corners(numberOfCorrespondences, 4);
    for( int i (0); i < numberOfCorrespondences; i++)
    {
        ISProcessor2::Correspondence *tempCorrespondence = &correspondences[i];
        corners(i, 0) = tempCorrespondence->cornerA->x;
        corners(i, 1) = tempCorrespondence->cornerA->y;
        corners(i, 2) = tempCorrespondence->cornerB->x;
        corners(i, 3) = tempCorrespondence->cornerB->y;
    }
    
    delete isprocessor;
    delete [] cornerImage1;
    delete [] cornerImage2;

    delete [] image1;
    delete [] image2;

    return corners;
}
