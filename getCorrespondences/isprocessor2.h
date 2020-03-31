/****************************************************************************
**
** Copyright (C) 2013 ISL.
** All rights reserved.
**
** Contact: Gwenael Schmitt (ISL-AVP) <gwenael.schmitt@isl.eu>
**
** version:
** - 20131024 - Colin Fischer
** - 20160418 - Martin Rebert
**
****************************************************************************/

#ifndef ISProcessor2_H
#define ISProcessor2_H

#include "gsmean.h"
#include "cdcorner.h"
#include "integralimage.h"
#include "math.h"
#include "cstdlib"
#include "cfloat"
#include "smmintrin.h"
#include <cstring>
#include "iostream"
//#include <qmath.h>

//------------------------------------------------------------------------------
class CoreProcess
{
public:
    CoreProcess(int width, int height, int margin=1);
    const unsigned char* output() const { return m_output; }
    void setInput(const unsigned char *input=NULL) { m_input = input;}
    void setOutput(unsigned char *output) { m_output = output; }
    virtual void process() = 0;

protected:
    int m_width;
    int m_height;

    int m_margin;
    int m_input_offset;
    int m_row_size;
    int m_row_offset;
    int m_rows_number;

    const unsigned char *m_input;
    unsigned char *m_output;
};

//------------------------------------------------------------------------------
class HarrisProcess : public CoreProcess
{
public:
    HarrisProcess(int width, int height, int margin = 2);
    ~HarrisProcess();
    void setSkyHeight(int skyHeight);
    void setX2Input(float* x2Input) { m_x2Input = x2Input; }
    void setY2Input(float* y2Input) { m_y2Input = y2Input; }
    void setXYInput(float* xyInput) { m_xyInput = xyInput; }
    void setCornersBuffer(CDCornerData* cornersBuffer, int cornersNumber);
    void process();
    const float* harrisValues() const { return m_harrisValues; }
    //float harrisMax() const { return m_harrisMax; }
private:
    float* m_x2Input;
    float* m_y2Input;
    float* m_xyInput;
    float* m_cornersMax;

    float m_zoneWidth;
    float m_zoneHeight;

    float m_cornerCellWidth;
    float m_cornerCellHeight;
    int   m_skyHeight;

    int* m_cornersIndexes;
    int m_cornersNumber;
    int m_sqrtCornersNumber;
    int m_sqrtCornersNumberD2;
    int m_sqrtCornersNumberD2_2;

    CDCornerData* m_cornersBuffer;

    float* m_harrisValues;
    //float m_harrisMax;
};

//------------------------------------------------------------------------------
class SobelProcess : public CoreProcess
{
public:
    SobelProcess(int width, int height);
    ~SobelProcess();
    float* x2Output() const { return m_x2Output; }
    float* y2Output() const { return m_y2Output; }
    float* xyOutput() const { return m_xyOutput; }
    void process();
private:
    float* m_x2Output;
    float* m_y2Output;
    float* m_xyOutput;
};

//------------------------------------------------------------------------------
class ISProcessor2
{
public:
    struct Correspondence
    {
        CDCornerData* cornerA;
        CDCornerData* cornerB;
        int match;
        int status;
        bool fbCorner;
    };

    struct CCIndex
    {
        int match;
        int indexA;
        int indexB;
    };

public:
    ISProcessor2(int width, int height, int cornersNumber);
    ~ISProcessor2();
    void setCornersImage1(CDCornerData *corners);
    void setCornersImage2(CDCornerData *corners);
    void setImage1(const unsigned char *image, CDCornerData *cornersBuffer);
    void setImage2(const unsigned char *image, CDCornerData *cornersBuffer);

    int process();
    int cornersNumber() const {return m_cornersNumber;}
    int meanPixelDist() const {return m_meanPixelDist;}
    void findCorners(const unsigned char *image, CDCornerData *cornersBuffer, const unsigned char *meanImage);

    SobelProcess& sobelProcess() { return m_sobelProcess; }
    HarrisProcess& harrisProcess() { return m_harrisProcess; }

    // Stats methods.
    int consistentCornerSetNumber() const { return m_consistentCornerSetNumber; }
    int inconsistentCornerSetNumber() const { return m_inconsistentCornerSetNumber; }
    int geometricalCornerDiff() const { return m_geometricalCornerDiff; }
    int numberOfFBCorners() const { return m_numberOfFBCorners; }
    int numberOfFBCornersSelectedByGeometric() const { return m_numberOfFBCornersSelectedByGeometric; }
    inline int getNumberOfCorrespondences() const {return m_correspondencesNumber;}
    Correspondence* getCorrespondences() {return m_correspondences;}


    int skyHeight() const { return m_skyHeight; }
    void setSkyHeight(int height);

private:
    unsigned int cornerMatch(const CDCornerData &cornerA, const CDCornerData &cornerB, int /*distance*/) const;
    void createCornerDesc(CDCornerData &corner, const unsigned char *meanImage, int width, unsigned int *briefIndices) const;
    void findCorrespondences(CDCornerData *cornersA, CDCornerData* cornersB, Correspondence* correspondences, int number);
    int findCorrespondences2(CDCornerData *cornersA, CDCornerData* cornersB, Correspondence* correspondences, int number);
    static int compareCorrespondance(const void * a, const void * b);
    void computeMeanDist(Correspondence **correspondences, int cornersNumber);

    /** Check the consistance of the corners of three juxtaposed zones
      * @return return true if the corner are consistent.
      */
    bool checkConsistentCorners(Correspondence *c1, Correspondence *c2, Correspondence *c3)
    {
        // The method used is to check the sign of the third value of the cross product between two vectors u and v formed with the three corners:
        // sign of u1*v2 - u2*v1
        // This operation is make for the given corners and their correspondance corner.

        bool a = (c2->cornerA->x - c1->cornerA->x) * (c3->cornerA->y - c1->cornerA->y) - (c2->cornerA->y - c1->cornerA->y) * (c3->cornerA->x - c1->cornerA->x) > 0;
        bool b = (c2->cornerB->x - c1->cornerB->x) * (c3->cornerB->y - c1->cornerB->y) - (c2->cornerB->y - c1->cornerB->y) * (c3->cornerB->x - c1->cornerB->x) > 0;
        return (a && b) || (!a && !b);
    }

private:
    int m_imageWidth;
    int m_imageHeight;
    int m_skyHeight;

    //BRIEF
    unsigned int* m_briefIndices;

    int m_rand_seed;

    CDCornerData* m_cornersImage1;
    CDCornerData* m_cornersImage2;

    const unsigned char *m_image1;
    const unsigned char *m_image2;

    const unsigned char *m_meanImage1;
    const unsigned char *m_meanImage2;

    CDCornerData* m_cornersImage1Buffer;
    CDCornerData* m_cornersImage2Buffer;

    unsigned char *m_gradientImage;

    int m_cornersNumber;

    SobelProcess m_sobelProcess;
    HarrisProcess m_harrisProcess;

    int m_subCornersNumber, m_correspondencesNumber;

    float** m_dissemblance;

    CDCornerData** m_cornerMinArray;
    int *m_cornerFlag1;
    int *m_cornerFlag2;
    int *m_zonesIndexesArray;
    CCIndex* m_ccIndexesArray;
    Correspondence* m_correspondences;

    int m_meanPixelDist;

    GSFloatingMean<int> m_meanCorrespondencesDistance;

    // Stats member attributes.
    int m_consistentCornerSetNumber;
    int m_inconsistentCornerSetNumber;
    int m_geometricalCornerDiff;
    int m_numberOfFBCorners;
    int m_numberOfFBCornersSelectedByGeometric;

};

#endif // ISProcessor2_H

