/****************************************************************************
**
** Copyright (C) 2013 ISL.
** All rights reserved.
**
** Contact: Gwenael Schmitt (ISL-AVP) <gwenael.schmitt@isl.eu>
**
****************************************************************************/

#ifndef GSFloatingMean_H
#define GSFloatingMean_H

//------------------------------------------------------------------------------
template <class T> class GSFloatingMean
{
public:
    GSFloatingMean(int size, T value=0);
    GSFloatingMean(const GSFloatingMean<T> &m);
    ~GSFloatingMean();
    GSFloatingMean<T>& operator+=(const T &value);
    GSFloatingMean<T>& operator=(const T &value);
    GSFloatingMean<T>& operator=(const GSFloatingMean<T> &m);
    inline operator T();
    inline operator T() const;
    inline T value() const;
    void setValue(const T &value);
    void updateValue(const T &value);
    inline int size() const { return m_fullSize; }
private:
    T *m_buffer;
    int m_index;
    T m_sum;
    int m_fullSize;
    int m_size;
};

//------------------------------------------------------------------------------
template <typename T> class GSMinMaxFloatingMean
{
public:
    GSMinMaxFloatingMean(int size, T value=0);
    GSMinMaxFloatingMean(const GSMinMaxFloatingMean<T> &m);
    ~GSMinMaxFloatingMean();
    GSMinMaxFloatingMean<T>& operator+=(const T &value);
    GSMinMaxFloatingMean<T>& operator=(const T &value);
    GSMinMaxFloatingMean<T>& operator=(const GSMinMaxFloatingMean<T> &m);
    inline operator T() { return m_sum / m_size; }
    inline operator T() const { return m_size!=0 ? m_sum / m_size : m_sum; }
    inline T value() const { return m_size!=0 ? m_sum / m_size : m_sum; }
    inline T min() const { return m_min; }
    inline T max() const { return m_max; }
    inline T last() const { return m_last; }
    inline T var() const {
        if (m_size==0) return m_sum2-(m_sum*m_sum);
        T mean=m_sum/m_size;
        return (m_sum2/m_size) - (mean*mean);
    }

    void setValue(const T &value);
    inline void reset(){ m_size=0; }
    void updateValue(const T &value);
    inline int size() const { return m_fullSize; }

private:
    T *m_buffer;
    int m_index;
    T m_sum;
    int m_fullSize;
    int m_size;
    T m_sum2;
    T m_min;
    T m_max;
    T m_last;
};

//------------------------------------------------------------------------------
template <typename T> class GSMean
{
public:
    inline GSMean(T value=0) : m_sum(value), m_size(0) { }
    GSMean<T>& operator+=(const T &value);
    GSMean<T>& operator=(const T &value);
    inline operator T() { return m_size!=0 ? m_sum / m_size : m_sum; }
    inline operator T() const { return m_size!=0 ? m_sum / m_size : m_sum; }
    inline T value() const { return m_sum / m_size; }
    void setValue(const T &value);
    void updateValue(const T &value);
    inline int size() { return m_size; }
private:
    T m_sum;
    int m_size;
};
//------------------------------------------------------------------------------
template <typename T> class GSMinMaxMean
{
public:
    inline GSMinMaxMean(T value=0) : m_sum(value), m_min(value), m_max(value), m_size(0) { }
    GSMinMaxMean<T>& operator+=(const T &value);
    GSMinMaxMean<T>& operator=(const T &value);
    inline operator T() { return m_sum / m_size; }
    inline operator T() const { return m_size!=0 ? m_sum / m_size : m_sum; }
    inline T value() const { return m_size!=0 ? m_sum / m_size : m_sum; }
    inline T min() const { return m_min; }
    inline T max() const { return m_max; }
    inline T last() const { return m_last; }
    inline T var() const {
        if (m_size==0) return m_sum2-(m_sum*m_sum);
        T mean=m_sum/m_size;
        return (m_sum2/m_size) - (mean*mean);
    }

    void setValue(const T &value);
    inline void reset(){ m_size=0; }
    void updateValue(const T &value);
    inline int size() const { return m_size; }

private:
    T m_sum;
    T m_sum2;
    T m_min;
    T m_max;
    T m_last;
    int m_size;
};


//---- GSFloatingMean ----------------------------------------------------------
template <class T> GSFloatingMean<T>::GSFloatingMean(int size, T value)
{
    m_fullSize = size > 0 ? size : 1;
    m_size = 0;
    m_buffer = new T[m_fullSize];
    m_index = 0;
    m_sum = value;
}

template <class T> GSFloatingMean<T>::GSFloatingMean(const GSFloatingMean<T> &m)
{
    m_fullSize = m.m_fullSize;
    m_size = m.m_size;
    m_buffer = new T[m_fullSize];
    m_index = m.m_index;
    m_sum = m.m_sum;
    memcpy(m_buffer, m.m_buffer, sizeof(T) * m_fullSize);
}

template <class T> GSFloatingMean<T>::~GSFloatingMean()
{
    delete[] m_buffer;
}

template <class T> GSFloatingMean<T>& GSFloatingMean<T>::operator+=(const T &value)
{
    if (m_size < m_fullSize)
    {
        if (m_size==0)
            m_sum = value;
        else
            m_sum += value;
        m_buffer[m_index] = value;
        m_size++;
    }
    else
    {
        m_sum += value;
        m_sum -= m_buffer[m_index];
        m_buffer[m_index] = value;
    }

    m_index++;
    if (m_index == m_fullSize)
        m_index = 0;

    return *this;
}

template <class T> GSFloatingMean<T>& GSFloatingMean<T>::operator=(const T &value)
{
    m_size = m_fullSize;
    m_sum = value * m_size;
    for (int i=0; i<m_size; i++)
        m_buffer[i] = value;
    return *this;
}

template <class T> GSFloatingMean<T>& GSFloatingMean<T>::operator=(const GSFloatingMean<T> &m)
{
    m_fullSize = m.m_fullSize;
    m_size = m.m_size;
    delete[] m_buffer;
    m_buffer = new T[m_fullSize];
    m_index = m.m_index;
    m_sum = m.m_sum;
    memcpy(m_buffer, m.m_buffer, sizeof(T) * m_fullSize);
    return *this;
}

template <class T> void GSFloatingMean<T>::setValue(const T &value)
{
    m_size = m_fullSize;
    m_sum = value * m_size;
    for (int i=0; i<m_size; i++)
        m_buffer[i] = value;
}

template <class T> void GSFloatingMean<T>::updateValue(const T &value)
{
    if (m_size < m_fullSize)
    {
        if (m_size==0)
            m_sum = value;
        else
            m_sum += value;
        m_buffer[m_index] = value;
        m_size++;
    }
    else
    {
        m_sum += value;
        m_sum -= m_buffer[m_index];
        m_buffer[m_index] = value;
    }

    m_index++;
    if (m_index == m_fullSize)
        m_index = 0;
}

template <class T> GSFloatingMean<T>::operator T()
{
    if (m_size==0)
        return m_sum;
    return m_sum / m_size;
}

template <class T> GSFloatingMean<T>::operator T() const
{
    if (m_size==0)
        return m_sum;
    return m_sum / m_size;
}

template <class T> T GSFloatingMean<T>::value() const
{
    if (m_size==0)
        return m_sum;
    return m_sum / m_size;
}

//---- GSMinMaxFloatingMean ----------------------------------------------------
template <class T> GSMinMaxFloatingMean<T>::GSMinMaxFloatingMean(int size, T value)
{
    m_fullSize = size > 0 ? size : 1;
    m_size = 0;
    m_buffer = new T[m_fullSize];
    m_index = 0;
    m_sum = value;
    m_sum2 = value*value;
    m_min = value;
    m_max = value;
}

template <class T> GSMinMaxFloatingMean<T>::GSMinMaxFloatingMean(const GSMinMaxFloatingMean<T> &m)
{
    m_fullSize = m.m_fullSize;
    m_size = m.m_size;
    m_buffer = new T[m_fullSize];
    m_index = m.m_index;
    m_sum = m.m_sum;
    m_sum2 = m.m_sum2;
    m_min = m.m_min;
    m_max =  m.m_max;
    memcpy(m_buffer, m.m_buffer, sizeof(T) * m_fullSize);
}

template <class T> GSMinMaxFloatingMean<T>::~GSMinMaxFloatingMean()
{
    delete[] m_buffer;
}

template <class T> GSMinMaxFloatingMean<T>& GSMinMaxFloatingMean<T>::operator+=(const T &value)
{
    if (m_size < m_fullSize)
    {
        if (m_size==0)
        {
            m_sum = value;
            m_sum2 = value*value;
        }
        else
        {
            m_sum += value;
            m_sum2 += value*value;
        }
        m_buffer[m_index] = value;
        m_size++;
    }
    else
    {
        m_sum += value;
        m_sum2 += value*value;
        m_sum -= m_buffer[m_index];
        m_buffer[m_index] = value;
    }
    if (value < m_min)
        m_min = value;
    if (value > m_max)
        m_max = value;

    m_index++;
    if (m_index == m_fullSize)
        m_index = 0;

    return *this;
}

template <class T> GSMinMaxFloatingMean<T>& GSMinMaxFloatingMean<T>::operator=(const T &value)
{
    m_size = m_fullSize;
    m_sum = value * m_size;
    m_sum2 = value*value*m_size;
    m_min= value;
    m_max = value;
    for (int i=0; i<m_size; i++)
        m_buffer[i] = value;
    return *this;
}

template <class T> GSMinMaxFloatingMean<T>& GSMinMaxFloatingMean<T>::operator=(const GSMinMaxFloatingMean<T> &m)
{
    m_fullSize = m.m_fullSize;
    m_size = m.m_size;
    delete[] m_buffer;
    m_buffer = new T[m_fullSize];
    m_index = m.m_index;
    m_sum = m.m_sum;
    m_sum2 = m.m_sum2;
    m_max = m.m_max;
    m_min = m.m_min;
    memcpy(m_buffer, m.m_buffer, sizeof(T) * m_fullSize);
    return *this;
}

template <class T> void GSMinMaxFloatingMean<T>::setValue(const T &value)
{
    m_size = m_fullSize;
    m_sum = value * m_size;
    m_sum2 = value*value*m_sum2;
    m_min = value;
    m_max = value;
    for (int i=0; i<m_size; i++)
        m_buffer[i] = value;
}

template <class T> void GSMinMaxFloatingMean<T>::updateValue(const T &value)
{
    if (m_size < m_fullSize)
    {
        if (m_size==0)
        {
            m_sum = value;
            m_sum2 = value*value;
        }
        else
        {
            m_sum += value;
            m_sum2 += value*value;
        }
        m_buffer[m_index] = value;
        m_size++;
    }
    else
    {
        m_sum += value;
        m_sum2 += value*value;
        m_sum -= m_buffer[m_index];
        m_buffer[m_index] = value;
    }
    if (value < m_min)
        m_min = value;
    if (value > m_max)
        m_max = value;

    m_index++;
    if (m_index == m_fullSize)
        m_index = 0;
}


//---- GSMean ------------------------------------------------------------------
template <typename T>
inline GSMean<T>& GSMean<T>::operator+=(const T &value)
{
    if (m_size==0)
        m_sum = value;
    else
        m_sum += value;
    m_size++;
    return *this;
}

template <typename T>
inline GSMean<T>& GSMean<T>::operator=(const T &value)
{
    m_sum = value;
    m_size = 1;
    return *this;
}

template <typename T>
inline void GSMean<T>::setValue(const T &value)
{
    m_sum = value;
    m_size = 1;
}

template <typename T>
inline void GSMean<T>::updateValue(const T &value)
{
    if (m_size==0)
        m_sum = value;
    else
        m_sum += value;
    m_size++;
}

//---- GSMinMaxMean ------------------------------------------------------------------
template <typename T>
inline GSMinMaxMean<T>& GSMinMaxMean<T>::operator+=(const T &value)
{
    if (m_size==0)
    {
        m_sum = value;
        m_sum2 = value*value;
        m_min = value;
        m_max = value;
    }
    else
    {
        m_sum += value;
        m_sum2 += value*value;
        if (value < m_min)
            m_min = value;
        if (value > m_max)
            m_max = value;
    }
    m_size++;
    m_last = value;
    return *this;
}

template <typename T>
inline GSMinMaxMean<T>& GSMinMaxMean<T>::operator=(const T &value)
{
    m_sum = value;
    m_sum2 = value*value;
    m_size = 1;
    m_min = value;
    m_max = value;
    m_last = value;
    return *this;
}

template <typename T>
inline void GSMinMaxMean<T>::setValue(const T &value)
{
    m_sum = value;
    m_sum2 = value*value;
    m_size = 1;
    m_min = value;
    m_max = value;
    m_last = value;
    return *this;
}

template <typename T>
inline void GSMinMaxMean<T>::updateValue(const T &value)
{
    if (m_size==0)
    {
        m_sum = value;
        m_sum2 = value*value;
        m_min = value;
        m_max = value;
    }
    else
    {
        m_sum += value;
        m_sum2 = value*value;
        if (value < m_min)
            m_min = value;
        if (value > m_max)
            m_max = value;
    }
    m_size++;
    m_last = value;
    return *this;
}

#endif // GSFloatingMean_H
