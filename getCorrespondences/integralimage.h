#ifndef INTEGRALIMAGE_H
#define INTEGRALIMAGE_H

class IntegralImage
{
public:
    IntegralImage(int width, int height);
    ~IntegralImage();
    void setImage(const unsigned char* imageBuffer);
    void meanFilter(int size);
    unsigned char* imageBuffer() const  { return m_imageBuffer; }
    unsigned int* integralImageBuffer() const { return m_integralImageBuffer; }

private:
    void createIntegralImage(const unsigned char* imageBuffer);
    void createMeanImage(const unsigned char *imageBuffer);

private:
    int m_width;
    int m_height;
    int m_size;
    int m_integralWidth;
    unsigned char* m_imageBuffer;
    unsigned int* m_integralImageBuffer;
};

#endif // INTEGRALIMAGE_H
