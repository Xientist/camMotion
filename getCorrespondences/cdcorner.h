#ifndef CDCORNER_H
#define CDCORNER_H

// BRIEF
#define CDCORNER_ROIWIDTH 64
#define CDCORNER_DESCSIZE 16
//static const int descCellWidth = 1;

struct CDCornerData
{
    float x;
    float y;
    unsigned int desc[CDCORNER_DESCSIZE];

    unsigned short refX;
    unsigned short refY;
    bool selected;
    bool fbCorner;
};

struct CDCornerData_d
{
    double x;
    double y;
//    quint32 desc[CDCORNER_DESCSIZE];

//    quint16 refX;
//    quint16 refY;
//    bool selected;
//    bool fbCorner;
};


#endif // CDCORNER_H
