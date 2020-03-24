#include <QCoreApplication>
#include <opencv2/opencv.hpp>
#include <string.h>

using namespace cv;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // construction des chemins vers les répertoires
    int indexSequence = 0;
    String s_index = (indexSequence>9)
            ? std::to_string(indexSequence)
            : "0"+std::to_string(indexSequence);

    String imageFolder = "dataset/sequences/" + s_index + "/image0/";
    String groundTruthFile = "poses/" + s_index + ".txt";

    // chargement des vérités terrain
    // à faire

    return a.exec();
}
