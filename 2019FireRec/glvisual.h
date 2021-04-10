#ifndef GLVISUAL_H
#define GLVISUAL_H

#include <QGLWidget>
#include <QVector2D>
#include <QtOpenGL>

extern std::ofstream stdof_LUTRGB;

const float pi = 3.141592653f;
const float twoPi = 2.0f * pi;
const float piBy2 = 0.5f * pi;
const float degToRad = pi / 180.0f;
const float radToDeg = 180.0f / pi;

class GLVisual : public QGLWidget {
Q_OBJECT

public:
    GLVisual( QWidget* parent = 0 );
    double rotate_y;
    double rotate_x;

protected:
    virtual void initializeGL();
    virtual void resizeGL( int w, int h );
    virtual void paintGL();

    virtual void keyPressEvent( QKeyEvent* e );
    virtual void mousePressEvent(QMouseEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);

private:
    float m_theta; /* < Rotation about x-axis */
    float m_phi;   /* < Rotation about y-axis */
    float m_aspectRatio;
    QVector2D m_lineWidthRange;
    float m_lineWidthStep;
    float m_lineWidth;
    QPoint lastPos;//11.3

};


#endif // GLVISUAL_H
