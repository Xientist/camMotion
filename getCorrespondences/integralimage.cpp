#include "integralimage.h"

IntegralImage::IntegralImage(int width, int height)
{
    m_width = width;
    m_height = height;
    m_size = width * height;
    m_integralWidth = width +1;
    m_imageBuffer = new unsigned char[width * height];
    m_integralImageBuffer = new unsigned int[(width+1) * (height+1)];
}

IntegralImage::~IntegralImage()
{
    delete[] m_imageBuffer;
    delete[] m_integralImageBuffer;
}

void IntegralImage::setImage(const unsigned char *imageBuffer)
{
    //m_imageBuffer = imageBuffer;
    createIntegralImage(imageBuffer);
}

void IntegralImage::meanFilter(int width)
{
    if (width > 1)
    {
        width |= 0x00000001;
        int widthD2 = width / 2;
        int squareWidth = width*width;
        for (int y = widthD2, ymax = m_height - widthD2; y < ymax; y++)
        {
            for (int x = widthD2, xmax = m_width - widthD2; x < xmax; x++)
            {
                m_imageBuffer[x + m_width * y] = ( m_integralImageBuffer[(x+widthD2+1) + m_integralWidth * (y+widthD2+1)]
                                                 - m_integralImageBuffer[(x+widthD2+1) + m_integralWidth * (y-widthD2)]
                                                 - m_integralImageBuffer[(x-widthD2) + m_integralWidth * (y+widthD2+1)]
                                                 + m_integralImageBuffer[(x-widthD2) + m_integralWidth * (y-widthD2)] ) / squareWidth;
            }
        }
    }
}

void IntegralImage::createIntegralImage(const unsigned char* imageBuffer)
{
    const unsigned char* inputP = imageBuffer;
    unsigned int* prevLineP = m_integralImageBuffer;
    unsigned char* outputP = m_imageBuffer;

    for (int i=0; i<=m_width; i++)
    {
        *prevLineP = 0;
    }

    unsigned int* integralP = m_integralImageBuffer + m_width + 1;
    prevLineP = m_integralImageBuffer + 1;
    *integralP = 0;
    integralP++;
    unsigned int sum = 0;
    for (int i=0, x=0; i<m_size; i++, x++)
    {
        *outputP = *inputP;
        if (x==m_width)
        {
            *integralP = 0;
            integralP++;
            prevLineP++;
            sum = 0;
            x = 0;
        }
        sum += *inputP;
        *integralP = *prevLineP + sum;
        integralP++;
        prevLineP++;
        inputP++;
        outputP++;
    }
}
