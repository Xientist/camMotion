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

#include "isprocessor2.h"

#define CDCORNER_ROIWIDTHD2 CDCORNER_ROIWIDTH / 2
#define CDCORNER_ROISIZE CDCORNER_ROIWIDTH * CDCORNER_ROIWIDTH
#define CDCORNER_MAXMATCH CDCORNER_ROISIZE * 255 + 1
#define CDIMAGE_MEANWIDTH 7

static unsigned short wordbits[65536];

//------------------------------------------------------------------------------
/*
      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 0  |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 1  |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 2  |   |   |(1)| # | # | # | # | # | # | # | # | # | # | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 3  |   |   | # | * | * | * | * | * | * | * | * | * | * | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 4  |   |   | # | * | * | * | * | * | * | * | * | * | * | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 5  |   |   | # | * | * | * | * | * | * | * | * | * | * | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 6  |   |   | # | * | * | * | * | * | * | * | * | * | * | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 7  |   |   | # | * | * | * | * | * | * | * | * | * | * | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 8  |   |   | # | * | * | * | * | * | * | * | * | * | * | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 9  |   |   | # | # | # | # | # | # | # | # | # | # | # | # |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
10  |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
11  |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

    m_margin = 3
    m_width = 16

    +---+---+---+---+---+---+---+---+---+---+
    | * | * | * | * | * | * | * | * | * | * | = row
    +---+---+---+---+---+---+---+---+---+---+
    |<--------- m_row_size = 10 ----------->|

    (1) = m_input_offset = 34

    +---+---+---+---+---+---+---+
    | # |   |   |   |   | # | * |
    +---+---+---+---+---+---+---+
    |<---- m_row_offset = 7 --->|

  */

CoreProcess::CoreProcess(int width, int height, int margin) :
    m_width(width), m_height(height), m_margin(margin), m_input(NULL), m_output(NULL)
{
    m_margin = margin < 1 ? 1 : margin;

    m_input_offset = m_margin - 1 + (m_margin - 1) * m_width;

    m_row_size = m_width - 2 * m_margin;
    m_row_offset = m_width - m_row_size + 1;

    m_rows_number = m_height - 2 * m_margin;
}
//------------------------------------------------------------------------------
HarrisProcess::HarrisProcess(int width, int height, int margin) :
    CoreProcess(width, height, margin < 2 ? 2 : margin),
    m_skyHeight(0)
{
    m_cornersNumber = 0;
    m_cornersMax = NULL;
    m_cornersBuffer = NULL;

    m_cornersIndexes = NULL;

    // Compensate the offset of the sobel output: one row before and one pixel before.
    m_input_offset -= (m_width + 1);

    m_harrisValues = new float[width * height];
    //m_harrisMax = 0;
}

HarrisProcess::~HarrisProcess()
{
    delete[] m_cornersMax;
    delete[] m_harrisValues;
    delete[] m_cornersIndexes;
}

void HarrisProcess::setSkyHeight(int skyHeight)
{
    if (skyHeight < m_margin)
        skyHeight = 0;
    m_input_offset += (skyHeight - m_skyHeight) * m_width;
    m_rows_number -= skyHeight - m_skyHeight;
    m_skyHeight = skyHeight;
}

//------------------------------------------------------------------------------
/*
    -- Example of buffer structure with 16 corners --

    Four zones of corners in the image:
    +---+---+ +---+---+
    |*0 |*1 | |*2 |*3 |
    +---+---+ +---+---+
    |*4 |*5 | |*6 |*7 |
    +---+---+ +---+---+
    +---+---+ +---+---+
    |*8 |*9 | |*10|*11|
    +---+---+ +---+---+
    |*12|*13| |*14|*15|
    +---+---+ +---+---+

    corner indexes array:
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    | 0 | 1 | 4 | 5 | 2 | 3 | 6 | 7 | 8 | 9 | 12| 13| 10| 11| 14| 15|
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

    corner buffer values after the execution of process() method:
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |*0 |*1 |*4 |*5 |*2 |*3 |*6 |*7 |*8 |*9 |*12|*13|*10|*11|*14|*15|
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    ^               ^               ^               ^               ^
    |    zone 0     |    zone 1     |    zone 2     |    zone 3     |
 */

void HarrisProcess::setCornersBuffer(CDCornerData *cornersBuffer, int cornersNumber)
{
    //qDebug()<<cornersNumber;
    if (m_cornersNumber != cornersNumber)
    {
        m_cornersNumber = cornersNumber;
        m_sqrtCornersNumber = sqrt(cornersNumber);
        m_sqrtCornersNumberD2 = m_sqrtCornersNumber / 2;
        m_sqrtCornersNumberD2_2 = m_sqrtCornersNumberD2 * m_sqrtCornersNumberD2;
        m_zoneWidth = static_cast<float>(m_row_size) / 2;
        m_zoneHeight = static_cast<float>(m_rows_number) / 2;

        m_cornerCellWidth = static_cast<float>(m_row_size) / m_sqrtCornersNumber;
        m_cornerCellHeight = static_cast<float>(m_rows_number) / m_sqrtCornersNumber;

        if (m_cornersMax)
            delete[] m_cornersMax;
        m_cornersMax = new float[m_cornersNumber];
    }

    m_cornersBuffer = cornersBuffer;

    // Clear old max value.
    for (int i=0; i<m_cornersNumber; i++)
    {
        m_cornersMax[i] = -FLT_MAX;
        std::memset(&m_cornersBuffer[i],0,sizeof(CDCornerData));
    }
}

void HarrisProcess::process()
{
    float* x2_in = m_x2Input + m_input_offset;
    float* y2_in = m_y2Input + m_input_offset;
    float* xy_in = m_xyInput + m_input_offset;

    float gx2, gy2, gxy;

    //m_harrisMax = 0;

    /*
      Gaussian filter:
      convolve the image with the horizontal matrix [1 2 1]
      then convolve the result with the vertical martix [1 2 1]
    */

    // Convolve the two first rows with the horizontal matrix:
    for (int i=0, imax=2*m_row_size, x=0; i<imax; i++)
    {
        gx2 = *x2_in + ((*(x2_in + 1)) * 2)  + *(x2_in + 2); // [1 2 1] Horizontal
        gy2 = *y2_in + ((*(y2_in + 1)) * 2)  + *(y2_in + 2); // [1 2 1] Horizontal
        gxy = *xy_in + ((*(xy_in + 1)) * 2)  + *(xy_in + 2); // [1 2 1] Horizontal

        *x2_in = gx2;
        *y2_in = gy2;
        *xy_in = gxy;

        x++;
        if (x < m_row_size)
        {
            x2_in++;
            y2_in++;
            xy_in++;
        }
        else // Go to next row.
        {
            x = 0;
            x2_in += m_row_offset;
            y2_in += m_row_offset;
            xy_in += m_row_offset;
        }
    }

    // Convolve with horizontal matrix and vertical matrix:
    for (int i=0, imax=m_row_size*m_rows_number, x=0, y=0; i<imax; i++)
        //for (int i=0, imax=roi.width()*roi.height(), x=1, y=1; i<1280*4; i++, x++)
    {
        // Horizontal convolution:
        gx2 = *x2_in + ((*(x2_in + 1)) * 2)  + *(x2_in + 2); // [1 2 1] Horizontal
        gy2 = *y2_in + ((*(y2_in + 1)) * 2)  + *(y2_in + 2); // [1 2 1] Horizontal
        gxy = *xy_in + ((*(xy_in + 1)) * 2)  + *(xy_in + 2); // [1 2 1] Horizontal

        *x2_in = gx2;
        *y2_in = gy2;
        *xy_in = gxy;

        // Vertical convolution:
        gx2 = *(x2_in - m_width * 2) + *(x2_in - m_width) * 2 + gx2; // [ 1 2 1] Vertical
        gy2 = *(y2_in - m_width * 2) + *(y2_in - m_width) * 2 + gy2; // [ 1 2 1] Vertical
        gxy = *(xy_in - m_width * 2) + *(xy_in - m_width) * 2 + gxy; // [ 1 2 1] Vertical


        // Harris is processed as : (E_x)^2*(E_y)^2 - E_xy*E_xy - (0.2*((E_x)^2+(E_y)^2))^2
        // i.e. with using the present variable as : gx2*gy2 - gxy*gxy - (gx2+gy2)*(gx2+gy2)*0.04;

        float harrisValue = gx2 * gy2 - gxy * gxy - (gx2 + gy2) * (gx2 + gy2) * 0.04;
        //if (harrisValue<0) harrisValue = -harrisValue;
        // Negative values correspond to border detection
        // and positive values correspond to corner detection
        //if (harrisValue<0) harrisValue = 0.0;

        m_harrisValues[m_margin + x + (m_margin + m_skyHeight + y)*m_width] = harrisValue;

        x++;
        if (x < m_row_size)
        {
            x2_in++;
            y2_in++;
            xy_in++;
        }
        else // Go to next row.
        {
            x = 0;
            y++;
            x2_in += m_row_offset;
            y2_in += m_row_offset;
            xy_in += m_row_offset;
        }
    }

    // 0 1 2
    // 7 V 3
    // 6 5 4

    float* harrisV = m_harrisValues + m_margin + (m_margin + m_skyHeight)*m_width;
    float* harris7 = harrisV - 1;
    float* harris3 = harrisV + 1;
    float* harris1 = harrisV - m_width;
    float* harris0 = harris1 - 1;
    float* harris2 = harris1 + 1;
    float* harris5 = harrisV + m_width;
    float* harris6 = harris5 - 1;
    float* harris4 = harris5 + 1;

    float sum, deltaX, deltaY;

    for (int i=0, imax=m_row_size*m_rows_number, x=0, y=0; i<imax; i++)
    {
        // Is local maximum?
        if (*harrisV > *harris0 && *harrisV > *harris1 && *harrisV > *harris2 && *harrisV > *harris3 && *harrisV > *harris4 && *harrisV > *harris5 && *harrisV > *harris6 && *harrisV > *harris7)
        {
//            if (*harrisV > m_harrisMax)
//                m_harrisMax = *harrisV;

            int cornerID = static_cast<int>(x / m_cornerCellWidth) + m_sqrtCornersNumber * static_cast<int>(y / m_cornerCellHeight);
            if (*harrisV >= m_cornersMax[cornerID])
            {
                sum = 0;
                (*harrisV > 0) ? sum += *harrisV : sum += - *harrisV;
                (*harris0 > 0) ? sum += *harris0 : sum += - *harris0;
                (*harris1 > 0) ? sum += *harris1 : sum += - *harris1;
                (*harris2 > 0) ? sum += *harris2 : sum += - *harris2;
                (*harris3 > 0) ? sum += *harris3 : sum += - *harris3;
                (*harris4 > 0) ? sum += *harris4 : sum += - *harris4;
                (*harris5 > 0) ? sum += *harris5 : sum += - *harris5;
                (*harris6 > 0) ? sum += *harris6 : sum += - *harris6;
                (*harris7 > 0) ? sum += *harris7 : sum += - *harris7;
                // 0 coordinates (-1,-1)
                // 1 coordinates (0,-1)
                // 2 coordinates (1,-1)
                // 3 coordinates (1,0)
                // 4 coordinates (1,1)
                // 5 coordinates (0,1)
                // 6 coordinates (-1,1)
                // 7 coordinates (-1,0)
                deltaX = *harris2 + *harris3 + *harris4 - *harris0 - *harris6 - *harris7;
                deltaX /= sum;
                deltaY = *harris4 + *harris5 + *harris6 - *harris0 - *harris1 - *harris2;
                deltaY /= sum;
                // Tej les deltas pour avoir du pixelll
                m_cornersMax[cornerID] = *harrisV;
                m_cornersBuffer[cornerID].x = m_margin + x + deltaX;
                m_cornersBuffer[cornerID].y = m_margin + m_skyHeight + y + deltaY;
            }
        }
        else
        {
            //m_input[m_margin + x + (m_margin + m_skyHeight + y)*m_width] = 0;
        }

        // Go to next value
        x++;
        if (x < m_row_size)
        {
            harrisV++;
            harris0++;
            harris1++;
            harris2++;
            harris3++;
            harris4++;
            harris5++;
            harris6++;
            harris7++;
        }
        else // Go to next row.
        {
            x = 0;
            y++;

            harrisV += m_row_offset;
            harris0 += m_row_offset;
            harris1 += m_row_offset;
            harris2 += m_row_offset;
            harris3 += m_row_offset;
            harris4 += m_row_offset;
            harris5 += m_row_offset;
            harris6 += m_row_offset;
            harris7 += m_row_offset;
        }
    }
}

//------------------------------------------------------------------------------
SobelProcess::SobelProcess(int width, int height) :
    CoreProcess(width, height)
{
    m_x2Output = new float[width*height];
    m_y2Output = new float[width*height];
    m_xyOutput = new float[width*height];

    for (int i=0, size = width*height; i<size; i++)
    {
        m_x2Output[i] = 0;
        m_y2Output[i] = 0;
        m_xyOutput[i] = 0;
    }
}

SobelProcess::~SobelProcess()
{
    delete[] m_x2Output;
    delete[] m_y2Output;
    delete[] m_xyOutput;
}

void SobelProcess::process()
{
    const unsigned char* in = m_input + m_input_offset;
    float* x_out = m_x2Output + m_input_offset;
    float* y_out = m_y2Output + m_input_offset;
    float* x2_out = x_out;
    float* y2_out = y_out;
    float* xy_out = m_xyOutput + m_input_offset;

    float gx, gy;

    /*
      Sobel x:
      convolve the image with the horizontal matrix [-1 0 1]
      then convolve the result with the vertical martix [1 2 1]

      Sobel y:
      convolve the image with the horizontal matrix [1 2 1]
      then convolve the result with the vertical martix [-1 0 1]
    */

    // Convolve the two first rows with the horizontal matrix:
    for (int i=0, imax=2*m_row_size, x=0; i<imax; i++)
    {
        gx = *(in + 2) - (*in);                    // [-1 0 1] Horizontal
        gy = *in + ((*(in + 1)) * 2)  + *(in + 2); // [ 1 2 1] Horizontal

        *x_out = gx;
        *y_out = gy;

        x++;
        if (x < m_row_size)
        {
            in++;
            x_out++;
            y_out++;
        }
        else // Go to next row.
        {
            x = 0;
            in += m_row_offset;
            x_out += m_row_offset;
            y_out += m_row_offset;
        }
    }

    // Convolve with horizontal matrix and vertical matrix:
    for (int i=0, imax=m_row_size*m_rows_number, x=0; i<imax; i++)
    {
        // Horizontal convolution:
        gx = *(in + 2) - (*in);                    // [-1 0 1] Horizontal
        gy = *in + ((*(in + 1)) * 2)  + *(in + 2); // [ 1 2 1] Horizontal

        *x_out = gx;
        *y_out = gy;


        // Vertical convolution:
        gx = *x2_out + *(x2_out + m_width) * 2 + gx; // [ 1 2 1] Vertical
        gy = gy - *y2_out;                           // [-1 0 1] Vertical

        *x2_out = gx * gx;
        *y2_out = gy * gy;
        *xy_out = gx * gy;

        x++;
        if (x < m_row_size)
        {
            in++;
            x_out++;
            y_out++;
            x2_out++;
            y2_out++;
            xy_out++;
        }
        else // Go to next row.
        {
            x = 0;
            in += m_row_offset;
            x_out += m_row_offset;
            y_out += m_row_offset;
            x2_out += m_row_offset;
            y2_out += m_row_offset;
            xy_out += m_row_offset;
        }
    }
}

//GSVector8 ISProcessor2__h, ISProcessor2__b;
//GSMatrix8x8 ISProcessor2__fastA;

//------------------------------------------------------------------------------
ISProcessor2::ISProcessor2(int width, int height, int cornersNumber):
    m_imageWidth(width), m_imageHeight(height),
    m_cornersNumber(cornersNumber),
    m_sobelProcess(width, height),
    m_harrisProcess(width, height, CDCORNER_ROIWIDTHD2),
    m_meanCorrespondencesDistance(10)
{
    m_subCornersNumber = m_cornersNumber/4;
    //m_subCornersNumber = m_cornersNumber;
    m_correspondencesNumber=m_subCornersNumber;

    m_dissemblance = new float*[m_subCornersNumber];
    for (int i=0; i<m_subCornersNumber; i++)
    {
        m_dissemblance[i] = new float[m_subCornersNumber];
    }

    m_skyHeight = 0;
    m_harrisProcess.setSkyHeight(m_skyHeight);
    m_harrisProcess.setX2Input(m_sobelProcess.x2Output());
    m_harrisProcess.setY2Input(m_sobelProcess.y2Output());
    m_harrisProcess.setXYInput(m_sobelProcess.xyOutput());

    m_cornerMinArray = new CDCornerData*[m_subCornersNumber];
    m_ccIndexesArray = new CCIndex[m_cornersNumber*m_cornersNumber];
    m_zonesIndexesArray = new int[m_cornersNumber];

    m_cornerFlag1 = new int[m_subCornersNumber];
    m_cornerFlag2 = new int[m_subCornersNumber];
    m_correspondences = new Correspondence[m_subCornersNumber];
    m_rand_seed = 123456;

    // Stats
    m_consistentCornerSetNumber = 0;
    m_inconsistentCornerSetNumber = 0;
    m_geometricalCornerDiff = 0;
    m_numberOfFBCorners = 0;
    m_numberOfFBCornersSelectedByGeometric = 0;

    //---- Create brief random sequence ----
    srand(m_rand_seed);

    int size = CDCORNER_DESCSIZE * 64; // descSize * sizeof unsigned int32 * 2
    m_briefIndices = new unsigned int[size];

    for (int i=0; i<size; i++)
    {
        bool flag = false;
        unsigned int index = 0;
        do
        {
            index = rand() % CDCORNER_ROIWIDTH + (rand() % CDCORNER_ROIWIDTH) * m_imageWidth;
            flag = false;
            for (int j=0; j<i; j++)
                if (m_briefIndices[j] == index)
                    flag = true;
        } while (flag);
        m_briefIndices[i] = index;
    }


    //---- Initialize bit count table if _mm_popcnt_u32 not used ----
    for (unsigned int i=0; i<65536; i++)
        wordbits[i] = __builtin_popcount(i);

    //---- Creation of the map index -> index begin zone
    int numberOfCornersByline=sqrt(m_cornersNumber);
    for(int i=0;i<m_cornersNumber;i++)
    {
        m_zonesIndexesArray[i]=i/(numberOfCornersByline*2)*numberOfCornersByline/2+(i%numberOfCornersByline)/2;
    }

    m_meanCorrespondencesDistance = 0;

    m_cornersImage1 = NULL;
    m_cornersImage2 = NULL;
    m_image1 = NULL;
    m_image2 = NULL;

}
ISProcessor2::~ISProcessor2()
{
    for (int i=0; i<m_subCornersNumber; i++)
    {
        delete[] m_dissemblance[i];
    }
    delete[] m_dissemblance;

    delete[] m_cornerMinArray;
    delete[] m_ccIndexesArray;
    delete[] m_zonesIndexesArray;
    delete[] m_cornerFlag1;
    delete[] m_cornerFlag2;
    delete[] m_correspondences;
    delete[] m_briefIndices;
}

void ISProcessor2::setSkyHeight(int height)
{
    m_skyHeight = height;
    m_harrisProcess.setSkyHeight(m_skyHeight);
}

void ISProcessor2::setCornersImage2(CDCornerData *corners)
{
    m_cornersImage2 = corners;
    m_image2 = NULL;
}

void ISProcessor2::setCornersImage1(CDCornerData *corners)
{
    m_cornersImage1 = corners;
    m_image1 = NULL;
}

void ISProcessor2::setImage2(const unsigned char *image, CDCornerData *cornersBuffer)
{
    m_image2 = image;
    m_cornersImage2 = NULL;
    m_cornersImage2Buffer = cornersBuffer;
}

void ISProcessor2::setImage1(const unsigned char *image, CDCornerData *cornersBuffer)
{
    m_image1 = image;
    m_cornersImage1 = NULL;
    m_cornersImage1Buffer = cornersBuffer;
}

void ISProcessor2::createCornerDesc(CDCornerData &corner, const unsigned char *meanImage, int width, unsigned int *briefIndices) const
{
    unsigned int *briefIndices2 = briefIndices + CDCORNER_DESCSIZE * 32;

    const unsigned char* imageIndex = meanImage + (unsigned short)corner.x - CDCORNER_ROIWIDTHD2 + ((unsigned short)corner.y - CDCORNER_ROIWIDTHD2) * width;

    unsigned int *d = corner.desc;

    for (int i=0; i<CDCORNER_DESCSIZE; i++)
    {
        *d = 0x00000000;

        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x80000000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x40000000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x20000000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x10000000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x08000000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x04000000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x02000000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x01000000;

        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00800000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00400000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00200000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00100000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00080000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00040000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00020000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00010000;

        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00008000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00004000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00002000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00001000;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000800;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000400;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000200;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000100;

        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000080;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000040;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000020;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000010;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000008;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000004;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000002;
        if (*(imageIndex + *briefIndices++) > *(imageIndex + *briefIndices2++)) *d |= 0x00000001;

        d++;
    }
}

unsigned int ISProcessor2::cornerMatch(const CDCornerData &cornerA, const CDCornerData &cornerB, int /*distance*/) const
{
    unsigned int result = 0;

    //---- BRIEF descriptor matching -----
    for (int i=0; i<CDCORNER_DESCSIZE; i++)
    {
        unsigned int v = cornerA.desc[i] ^ cornerB.desc[i];
        result += _mm_popcnt_u32(v);
    }

    //---- Distance weight ----
    /*int dx = m_x - corner.m_x;
    int dy = m_y - corner.m_y;
    result *= 1.0 + 0.00000390625 * (dx*dx + dy*dy);*/

    int dx = cornerA.x - cornerB.x;
    int dy = cornerA.y - cornerB.y;
    int d = dx*dx + dy*dy;
    //if ((d < distance-100) || (distance+100 < d))
    //if (d > distance+100)
    if (d>409600)
        result += 512;

    return result;
}

void ISProcessor2::findCorners(const unsigned char *image, CDCornerData *cornersBuffer, const unsigned char *meanImage)
{
    m_sobelProcess.setInput(image);
    m_sobelProcess.process();

    //m_harrisProcess.setInput(m_integralImage.imageBuffer());
    m_harrisProcess.setInput(image);
    m_harrisProcess.setCornersBuffer(cornersBuffer, m_cornersNumber);
    m_harrisProcess.process();
    //qDebug()<<"ISProcessor2::findCorners****************************new loop **************************";
    //for (int i=0; i<m_cornersNumber; i++)
    //    qDebug()<<"ISProcessor2::findCorners"<<cornersBuffer[i].x<<cornersBuffer[i].y;

    for (int i=0; i<m_cornersNumber; i++)
    {
        createCornerDesc(cornersBuffer[i], meanImage, m_imageWidth, m_briefIndices);
    }
}

int ISProcessor2::compareCorrespondance(const void * a, const void * b)
{
    return (static_cast<const CCIndex*>(a))->match - (static_cast<const CCIndex*>(b))->match;
}
//----FB selection
int ISProcessor2::findCorrespondences2(CDCornerData *cornersA, CDCornerData *cornersB, Correspondence *correspondences, int number)
{
    int numberOfCornersByline = sqrt(m_cornersNumber);
    CDCornerData *cornerA;
    CDCornerData *cornerB;
    float match;
    Correspondence *correspondence;

    CCIndex* pCCIndex = m_ccIndexesArray;
    for (int i=0; i<number; i++)
    {
        cornerA = &cornersA[i];
        {
            for (int j=0; j<number; j++)
            {
                cornerB = &cornersB[j];
                match = cornerMatch(*cornerA, *cornerB, m_meanCorrespondencesDistance.value());
                pCCIndex->match = match;
                pCCIndex->indexA = i;
                pCCIndex->indexB = j;
                pCCIndex++;

            }
        }
    }
    for(int i=0;i<m_subCornersNumber;i++)
    {
        m_cornerFlag1[i] = 0;
        m_cornerFlag2[i] = 0;
    }


    int indexCorrespondances = 0;

    int indexAMin1 = 0;
    int indexAMin2 = 0;
    int indexAMin3 = 0;
    int indexAMin4 = 0;
    for(int i=0,c=0;i<m_subCornersNumber;i++)
    {
        // search of the min of the row
        CCIndex *pCCIndex1 = m_ccIndexesArray + c;
        if(!m_cornerFlag2[m_zonesIndexesArray[pCCIndex1->indexB]])
        {
            CCIndex *pCCIndex2 = pCCIndex1+1;
            CCIndex *pCCIndex3 = pCCIndex1+numberOfCornersByline;
            CCIndex *pCCIndex4 = pCCIndex3+1;

            int indexB1 = pCCIndex1->indexB;
            int indexB2 = pCCIndex2->indexB;
            int indexB3 = pCCIndex3->indexB;
            int indexB4 = pCCIndex4->indexB;

            int min1 = CDCORNER_MAXMATCH;
            int min2 = CDCORNER_MAXMATCH;
            int min3 = CDCORNER_MAXMATCH;
            int min4 = CDCORNER_MAXMATCH;

            for(int j=0;j<m_cornersNumber;j++)
            {
                if(pCCIndex1->match<min1)
                {
                    min1 = pCCIndex1->match;
                    indexAMin1 = j;
                }
                if(pCCIndex2->match<min2)
                {
                    min2=pCCIndex2->match;
                    indexAMin2 = j;
                }
                if(pCCIndex3->match<min3)
                {
                    min3=pCCIndex3->match;
                    indexAMin3 = j;
                }
                if(pCCIndex4->match<min4)
                {
                    min4=pCCIndex4->match;
                    indexAMin4 = j;
                }
                pCCIndex1+=m_cornersNumber;
                pCCIndex2+=m_cornersNumber;
                pCCIndex3+=m_cornersNumber;
                pCCIndex4+=m_cornersNumber;
            }

            // see if it is also the min of the line
            pCCIndex1 = m_ccIndexesArray + indexAMin1 * m_cornersNumber;
            pCCIndex2 = m_ccIndexesArray + indexAMin2 * m_cornersNumber;
            pCCIndex3 = m_ccIndexesArray + indexAMin3 * m_cornersNumber;
            pCCIndex4 = m_ccIndexesArray + indexAMin4 * m_cornersNumber;
            bool trueMin1=true;
            bool trueMin2=true;
            bool trueMin3=true;
            bool trueMin4=true;
            for(int j=0;j<m_cornersNumber;j++)
            {
                if(pCCIndex1->match<min1) // the min of the line is not the min of the row
                {
                    j=m_cornersNumber; // we stop here and don't consider it
                    trueMin1=false;
                }
                pCCIndex1++;
            }
            for(int j=0;j<m_cornersNumber;j++)
            {
                if(pCCIndex2->match<min2) // the min of the line is not the min of the row
                {
                    j=m_cornersNumber; // we stop here and don't consider it
                    trueMin2=false;
                }
                pCCIndex2++;
            }
            for(int j=0;j<m_cornersNumber;j++)
            {
                if(pCCIndex3->match<min3) // the min of the line is not the min of the row
                {
                    j=m_cornersNumber; // we stop here and don't consider it
                    trueMin3=false;
                }
                pCCIndex3++;
            }
            for(int j=0;j<m_cornersNumber;j++)
            {
                if(pCCIndex4->match<min4) // the min of the line is not the min of the row
                {
                    j=m_cornersNumber; // we stop here and don't consider it
                    trueMin4=false;
                }
                pCCIndex4++;
            }

            bool trueMin = false;
            int min = CDCORNER_MAXMATCH;
            int indexA, indexB;
            if (trueMin1 && !m_cornerFlag1[m_zonesIndexesArray[indexAMin1]])
            {
                min = min1;
                trueMin = true;
                indexA = indexAMin1;
                indexB = indexB1;
            }
            if (trueMin2 && !m_cornerFlag1[m_zonesIndexesArray[indexAMin2]] && min2<min)
            {
                min = min2;
                trueMin = true;
                indexA = indexAMin2;
                indexB = indexB2;
            }
            if (trueMin3 && !m_cornerFlag1[m_zonesIndexesArray[indexAMin3]] && min3<min)
            {
                min = min3;
                trueMin = true;
                indexA = indexAMin3;
                indexB = indexB3;
            }
            if (trueMin4 && !m_cornerFlag1[m_zonesIndexesArray[indexAMin4]] && min4<min)
            {
                min = min4;
                trueMin = true;
                indexA = indexAMin4;
                indexB = indexB4;
            }

            if(trueMin)
            {

                // if it is a true min, and the indexes are not allready marked, we add it in correspondence
                //qDebug()<<"new correspondance " << indexCorrespondances << min;
                correspondence = &correspondences[indexCorrespondances];
                correspondence->cornerA = &cornersA[indexA];
                correspondence->cornerB = &cornersB[indexB];
                correspondence->match = min;
                correspondence->status = 0;
                correspondence->fbCorner = true;

                m_cornerFlag2[m_zonesIndexesArray[indexB]] = 1;
                m_cornerFlag1[m_zonesIndexesArray[indexA]] = 1;

                indexCorrespondances++;

            }
        }

        c+=2;
        if (c%numberOfCornersByline==0)
        {
            c+=numberOfCornersByline;
        }
    }
    return indexCorrespondances;
}

//---- Winner take all ----
void ISProcessor2::findCorrespondences(CDCornerData *cornersA, CDCornerData *cornersB, Correspondence *correspondences, int number)
{
    const CDCornerData *cornerA;
    CDCornerData *cornerB;
    int match;
    Correspondence *correspondence;

    CCIndex* pCCIndex = m_ccIndexesArray;

    for (int i=0; i<number; i++)
    {
        cornerA = &cornersA[i];
        if(cornerA->x!=CDCORNER_ROIWIDTHD2)
        {
            for (int j=0; j<number; j++)
            {
                //qDebug()<<"j = "<<j << match;
                cornerB = &cornersB[j];
                if(cornerB->x!=CDCORNER_ROIWIDTHD2)
                {
                    match = cornerMatch(*cornerA, *cornerB, m_meanCorrespondencesDistance.value());
                    m_cornerFlag2[i] = 0;
                }
                else
                {
                    match = CDCORNER_MAXMATCH;
                    m_cornerFlag2[i] = 0;
                }
                pCCIndex->match = match;
                pCCIndex->indexA = i;
                pCCIndex->indexB = j;
                pCCIndex++;

            }
            m_cornerFlag1[i] = 0;
        }
        else
            m_cornerFlag1[i] = 1;

        //qDebug()<<"i = "<<i << match;
        //m_cornerFlag1[i] = 0;
    }

    qsort(m_ccIndexesArray, m_cornersNumber*m_cornersNumber, sizeof(CCIndex), ISProcessor2::compareCorrespondance);

    pCCIndex = m_ccIndexesArray;
    int zoneIDA, zoneIDB;
    int numberOfCornersByline = sqrt(m_cornersNumber);
    // there are 1024 corners
    // 4 corners by zone and we should select only one corner by zone
    //
    for (int i=0; i<m_subCornersNumber; i++)
    {
        while(m_cornerFlag1[pCCIndex->indexA] || m_cornerFlag2[pCCIndex->indexB]) pCCIndex++;
        correspondence = &correspondences[i];
        zoneIDA = pCCIndex->indexA/(numberOfCornersByline*2)*numberOfCornersByline/2+(pCCIndex->indexA%numberOfCornersByline)/2; //
        zoneIDB = pCCIndex->indexB/(numberOfCornersByline*2)*numberOfCornersByline/2+(pCCIndex->indexB%numberOfCornersByline)/2; //
        //qDebug()<<"zone ID A : " <<zoneIDA;
        //qDebug()<<"zone ID B : " <<zoneIDB;
        correspondence->cornerA = &cornersA[pCCIndex->indexA];
        correspondence->cornerB = &cornersB[pCCIndex->indexB];
        correspondence->match = pCCIndex->match;
        correspondence->status = 0;
        // we have to marque all the zone, not only the corner
        int line = (2*zoneIDA)/numberOfCornersByline;
        int lineB = (2*zoneIDB)/numberOfCornersByline;
        m_cornerFlag1[2*zoneIDA+line*numberOfCornersByline] = 1;
        m_cornerFlag1[2*zoneIDA+line*numberOfCornersByline+1] = 1;
        m_cornerFlag1[2*zoneIDA+(line+1)*numberOfCornersByline] = 1;
        m_cornerFlag1[2*zoneIDA+(line+1)*numberOfCornersByline+1] = 1;

        m_cornerFlag2[2*zoneIDB+lineB*numberOfCornersByline] = 1;
        m_cornerFlag2[2*zoneIDB+lineB*numberOfCornersByline+1] = 1;
        m_cornerFlag2[2*zoneIDB+(lineB+1)*numberOfCornersByline] = 1;
        m_cornerFlag2[2*zoneIDB+(lineB+1)*numberOfCornersByline+1] = 1;
    }
}

void ISProcessor2::computeMeanDist(Correspondence **correspondences, int cornersNumber)
{
    int sum=0;
    for(int i=0;i<cornersNumber;i++)
    {
        sum+=correspondences[i]->cornerA->x-correspondences[i]->cornerB->x;
    }
    m_meanPixelDist=sum/cornersNumber;
}

int ISProcessor2::process()
{
    for (int i=0; i<m_cornersNumber; i++)
    {
        CDCornerData *cornerImage2 = &m_cornersImage2[i];
        //CDCornerData *cornerA = &m_cornersA[i];
        cornerImage2->refX = 0;
        cornerImage2->refY = 0;
        cornerImage2->selected = false;
    }

    m_correspondencesNumber = findCorrespondences2(m_cornersImage1 , m_cornersImage2, m_correspondences, m_cornersNumber);


    m_numberOfFBCorners=0;
    for (int i=0; i<m_correspondencesNumber; i++)
    {
        Correspondence *correspondence = &m_correspondences[i];
        correspondence->cornerB->refX = correspondence->cornerA->x;
        correspondence->cornerB->refY = correspondence->cornerA->y;
        if(correspondence->fbCorner)
        {
            m_numberOfFBCorners++;
            correspondence->cornerB->fbCorner = true;
        }
        else
            correspondence->cornerB->fbCorner = false;
    }

    return 0;
}

// ==========================================================================================================
