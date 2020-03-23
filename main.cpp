#include <QCoreApplication>
#include <opencv2/opencv.hpp>
#include <string.h>

using namespace cv;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    String imageName = "lena.png";

    std::cout << "lecture de l'image '" + imageName + "'...\n";
    Mat image = imread(imageName);

    namedWindow("test", WINDOW_AUTOSIZE);

    std::cout << "affichage de l'image...\n";
    imshow("test", image);
    std::cout << imageName + " affichée.\n";

    waitKey(0);

    destroyAllWindows();

    return a.exec();
}
