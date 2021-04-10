#include "glvisual.h"

#include <QCoreApplication>
#include <QDebug>
#include <QKeyEvent>
#include <math.h>
#include "mainwindow.h"

#include <QPainter>

std::ofstream stdof_LUTRGB;


GLVisual::GLVisual( QWidget* parent )
    : QGLWidget( parent ),
      m_theta( 0.0f ),
      m_phi( 0.0f ),
      m_aspectRatio( 1.0 ),
      m_lineWidthRange(),
      m_lineWidthStep( 0.0f ),
      m_lineWidth( 1.0f ) {
}


//初始化
void GLVisual::initializeGL() {
    // Set the clear color to black
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );

    // Set the drawing color to green
    //glColor3f( 0.0f, 1.0f, 0.0f );

    // Query some info about supported point sizes
    glGetFloatv( GL_LINE_WIDTH_RANGE, reinterpret_cast<float*>( &m_lineWidthRange ) );
    glGetFloatv( GL_LINE_WIDTH_GRANULARITY, &m_lineWidthStep );

    qDebug() << "Point size range:" << m_lineWidthRange;
    qDebug() << "Point size step:" << m_lineWidthStep;

    m_lineWidth = m_lineWidthRange.x();


    setFocusPolicy(Qt::StrongFocus);

}

void GLVisual::resizeGL( int w, int h ) {
    // Prevent a divde by zero
    if ( h == 0 )
        h = 1;

    // Set the viewport to window dimensions
    glViewport( 0, 0, w, h );

    // reset the coordinate system
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

    // Establish the clipping volume by setting up an orthographic projection
    double range = 100.0;
    m_aspectRatio = double( w ) / double( h );

    if ( w <= h )
        glOrtho( -range, range, -range / m_aspectRatio, range / m_aspectRatio, range, -range );
    else
        glOrtho( -range * m_aspectRatio, range * m_aspectRatio, -range, range, range, -range );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}

int openLUT = 0;
void GLVisual::paintGL() {
    // Clear the buffer with the current clearing color
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    if(0==openLUT){
        stdof_LUTRGB.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "LUT.csv", std::ios::out | std::ios::app);
        stdof_LUTRGB <<  "label\tiR\tiG\tiB\tpx\tpy\tpz\n";
        //stdof_LUTRGB.close();
        openLUT = 1;
    }

    // Set drawing colour to red (RGB)
    //glColor3f( 1.0f, 0.0f, 0.0f );

    // Save matrix state and do the custom rotation
    glPushMatrix();
    glRotatef( rotate_x, 1.0f, 0.0f, 0.0f );   //m_theta
    glRotatef( rotate_y, 0.0f, 1.0f, 0.0f ); //m_phi

    // Draw some Lines in a helix
    glLineWidth( m_lineWidth );
    glPointSize( m_lineWidth );

    //std::cout << "qs_fileName :" << qs_fileName.toStdString() << ")\n";
    /// Read
    std::ifstream LUTifR("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_R.txt");
    std::ifstream LUTifG("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_G.txt");
    std::ifstream LUTifB("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_B.txt");
    //std::ifstream LUTifLabel("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_Label.txt");

    std::string sR;
    std::string sG;
    std::string sB;
    //std::string sC;


    ///QFont font;
    ///QFont font("Serif", 16, QFont::Bold);


    //QFont font("Helvetica [Cronyx]", 13);
    //QFont font("San Serif", 18);
    ///font.setPointSize(14);

    qglColor(QColor("#0000ff"));
    //renderText(-50, -50, -50,"  Blue ", font);

    qglColor(QColor("#000000"));
    //renderText(-50, -50,  50,"  Black", font);

    qglColor(QColor("#00cccc"));
    //renderText(-50,  50, -50,"  Cyan", font);

    qglColor(QColor("#ff00ff"));
    //renderText( 50, -50, -50,"  Magenta", font);

    qglColor(QColor("#999999"));
    //renderText( 50,  50, -50,"  White", font);

    qglColor(QColor("#ff0000"));
    //renderText( 50, -50,  50,"  Red", font);

    qglColor(QColor("#00cc00"));
    //renderText(-50,  50,  50,"  Green", font);

    qglColor(QColor("#cccc00"));
    //renderText( 50,  50,  50,"  Yellow)", font);

    ///renderText(0, 0, 0,"(0, 0, 0)", font);

    glBegin(GL_POINTS);

        while ( std::getline( LUTifR, sR ) ){


            int R = std::stoi( sR.substr(0) , 0 );

            std::getline( LUTifG, sG );
            int G = std::stoi( sG.substr(0) , 0 );

            std::getline( LUTifB, sB );
            int B = std::stoi( sB.substr(0) , 0 );

            //printf("RGB components of the pixel selected: %d %d %d\n", R, G, B);

            //std::cout << "GLVisual: (R, G, B) = (" << sR << ", " << sG << ", " << sB << ")\n";

            glColor3f(R/256.0, G/256.0, B/256.0);
            glVertex3f( (R/256.0)*100.0-50, (G/256.0)*100.0-50.0, 50.0-(B/256.0)*100.0 );


            //std::getline( LUTifLabel, sC );
            //int label = std::stoi( sC.substr(0) , 0 );
            stdof_LUTRGB << "1\t" +
                            std::to_string( R ) + "\t" +
                            std::to_string( G ) + "\t" +
                            std::to_string( B ) + "\t" +
                            std::to_string( (R/256.0)*100.0-50 ) + "\t" +
                            std::to_string( (G/256.0)*100.0-50.0 ) + "\t" +
                            std::to_string( 50.0-(B/256.0)*100.0 ) + "\n";


        }

    glEnd();


    /// Face 1 Color RGB
    glBegin(GL_LINE_STRIP);
        /// #1 bk-g-y-r-bk
        glColor3f( 0.0, 0.0, 0.0 );       // P8 is blue --> BLACK
        glVertex3f( -50.0, -50.0,  50.0 );

        glColor3f( 0.0, 1.0, 0.0 );       // P7 is cyan --> GREEN ok
        glVertex3f( -50.0,  50.0,  50.0 );

        glColor3f( 1.0, 1.0, 0.0 );       // P6 is green --> YELLOW ok
        glVertex3f(  50.0,  50.0,  50.0 );

        glColor3f( 1.0, 0.0, 0.0 );       // P5 is black --> RED ok
        glVertex3f(  50.0, -50.0,  50.0 );

        glColor3f( 0.0, 0.0, 0.0 );       // P8 is blue --> BLACK ok
        glVertex3f( -50.0, -50.0,  50.0 );


        /// #2 bk-g-c-bl-bk
        //glColor3f( 0.0, 0.0, 0.0 );       // P8 is blue --> BLACK
        //glVertex3f( -50.0, -50.0,  50.0 );

        glColor3f( 0.0, 1.0, 0.0 );       // P7 is cyan --> GREEN
        glVertex3f( -50.0,  50.0,  50.0 );

        glColor3f( 0.0, 1.0, 1.0 );       // P3 is white --> CYAN ok
        glVertex3f( -50.0,  50.0, -50.0 );

        glColor3f( 0.0, 0.0, 1.0 );       // P4 is magenta --> BLUE ok
        glVertex3f( -50.0, -50.0, -50.0 );

        glColor3f( 0.0, 0.0, 0.0 );       // P8 is blue --> BLACK
        glVertex3f( -50.0, -50.0,  50.0 );


        /// #3 bk-r-y-w-m-r
        //glColor3f( 0.0, 0.0, 0.0 );       // P8 is blue --> BLACK
        //glVertex3f( -50.0, -50.0,  50.0 );

        glColor3f( 1.0, 0.0, 0.0 );       // P5 is black --> RED
        glVertex3f(  50.0, -50.0,  50.0 );

        glColor3f( 1.0, 1.0, 0.0 );       // P6 is green --> YELLOW
        glVertex3f(  50.0,  50.0,  50.0 );

        glColor3f( 1.0, 1.0, 1.0 );       // P2 is yellow --> WHITE ok
        glVertex3f(  50.0,  50.0, -50.0 );

        glColor3f( 1.0, 0.0, 1.0 );       // P1 is red --> MAGENTA ok
        glVertex3f(  50.0, -50.0, -50.0 );

        glColor3f( 1.0, 0.0, 0.0 );       // P5 is black --> RED
        glVertex3f(  50.0, -50.0,  50.0 );


        /// #4 r-y-g-c-w-y
        //glColor3f( 1.0, 0.0, 0.0 );       // P5 is black --> RED
        //glVertex3f(  50.0, -50.0,  50.0 );

        glColor3f( 1.0, 1.0, 0.0 );       // P6 is green --> YELLOW
        glVertex3f(  50.0,  50.0,  50.0 );

        glColor3f( 0.0, 1.0, 0.0 );       // P7 is cyan --> GREEN
        glVertex3f( -50.0,  50.0,  50.0 );

        glColor3f( 0.0, 1.0, 1.0 );       // P3 is white --> CYAN
        glVertex3f( -50.0,  50.0, -50.0 );

        glColor3f( 1.0, 1.0, 1.0 );       // P2 is yellow --> WHITE
        glVertex3f(  50.0,  50.0, -50.0 );

        glColor3f( 1.0, 1.0, 0.0 );       // P6 is green --> YELLOW
        glVertex3f(  50.0,  50.0,  50.0 );


        /// #5 y-r-m-bl-bk-r
        //glColor3f( 1.0, 1.0, 0.0 );       // P6 is green --> YELLOW
        //glVertex3f(  50.0,  50.0,  50.0 );

        glColor3f( 1.0, 0.0, 0.0 );       // P5 is black --> RED
        glVertex3f(  50.0, -50.0,  50.0 );

        glColor3f( 1.0, 0.0, 1.0 );       // P1 is red --> MAGENTA
        glVertex3f(  50.0, -50.0, -50.0 );

        glColor3f( 0.0, 0.0, 1.0 );       // P4 is magenta --> BLUE
        glVertex3f( -50.0, -50.0, -50.0 );

        glColor3f( 0.0, 0.0, 0.0 );       // P8 is blue --> BLACK
        glVertex3f( -50.0, -50.0,  50.0 );

        glColor3f( 1.0, 0.0, 0.0 );       // P5 is black --> RED
        glVertex3f(  50.0, -50.0,  50.0 );


        /// #6 r-m-bl-c-w-m

        /// #7 r-bk-w
        //glColor3f( 1.0, 0.0, 0.0 );       // P5 is black --> RED
        //glVertex3f(  50.0, -50.0,  50.0 );

        glColor3f( 0.0, 0.0, 0.0 );       // P8 is blue --> BLACK
        glVertex3f( -50.0, -50.0,  50.0 );

        //glColor3f( 1.0, 1.0, 1.0 );       // P2 is black --> WHITE
        //glVertex3f(  50.0,  50.0, -50.0 );


    glEnd();


/*
    // White Box
    glBegin(GL_LINE_STRIP);
        /// #1 bk-g-y-r-bk
        glColor3f( 0.7, 0.7, 0.7 );
        glVertex3f( -50.0, -50.0,  50.0 );
        glVertex3f( -50.0,  50.0,  50.0 );
        glVertex3f(  50.0,  50.0,  50.0 );
        glVertex3f(  50.0, -50.0,  50.0 );
        glVertex3f( -50.0, -50.0,  50.0 );

        glVertex3f( -50.0,  50.0,  50.0 );
        glVertex3f( -50.0,  50.0, -50.0 );
        glVertex3f( -50.0, -50.0, -50.0 );
        glVertex3f( -50.0, -50.0,  50.0 );

        glVertex3f(  50.0, -50.0,  50.0 );
        glVertex3f(  50.0,  50.0,  50.0 );
        glVertex3f(  50.0,  50.0, -50.0 );
        glVertex3f(  50.0, -50.0, -50.0 );
        glVertex3f(  50.0, -50.0,  50.0 );

        glVertex3f(  50.0,  50.0,  50.0 );
        glVertex3f( -50.0,  50.0,  50.0 );
        glVertex3f( -50.0,  50.0, -50.0 );
        glVertex3f(  50.0,  50.0, -50.0 );
        glVertex3f(  50.0,  50.0,  50.0 );

        glVertex3f(  50.0, -50.0,  50.0 );
        glVertex3f(  50.0, -50.0, -50.0 );
        glVertex3f( -50.0, -50.0, -50.0 );
        glVertex3f( -50.0, -50.0,  50.0 );
        glVertex3f(  50.0, -50.0,  50.0 );

        glVertex3f( -50.0, -50.0,  50.0 );
        //glVertex3f(  50.0,  50.0, -50.0 );

    glEnd();
*/


    glFlush();
    ///swapBuffers();

    // Restore the matrix state
    glPopMatrix();
    setAutoBufferSwap(true);

}

//键盘控制
void GLVisual::keyPressEvent( QKeyEvent* e ) {
    switch ( e->key() )
    {
    case Qt::Key_Escape:
        QCoreApplication::instance()->quit();
    break;

    case Qt::Key_Left:
        m_phi += 1.0f;
        rotate_y -= 5;
        updateGL();
    break;

    case Qt::Key_Right:
        m_phi -= 1.0f;
        rotate_y += 5;
        updateGL();
    break;

    case Qt::Key_Up:
        m_theta += 1.0f;
        rotate_x -= 5;
        updateGL();
    break;

    case Qt::Key_Down:
        m_theta -= 1.0f;
        rotate_x += 5;
        updateGL();
    break;

    case Qt::Key_Plus:
        m_lineWidth = qMin( m_lineWidth + m_lineWidthStep, float( m_lineWidthRange.y() ) );
        qDebug() << "m_lineWidth =" << m_lineWidth;
        updateGL();
    break;

    case Qt::Key_Minus:
        m_lineWidth = qMax( m_lineWidth - m_lineWidthStep, float( m_lineWidthRange.x() ) );
        qDebug() << "m_lineWidth =" << m_lineWidth;
        updateGL();
    break;

    default:
        //rotate_y -= 5;
        //rotate_x -= 5;
        updateGL();
        QGLWidget::keyPressEvent( e );
    }
}

//2016.11.3 21:17

void GLVisual::mousePressEvent(QMouseEvent *e){

    if (e->button()==Qt::LeftButton) {

        lastPos = e->pos();

        updateGL();
    }
}

void GLVisual::mouseMoveEvent(QMouseEvent *e){

    int dx = e->x() - lastPos.x();

    int dy = e->y() - lastPos.y();

    if (e->buttons()&Qt::LeftButton) {

        setMouseTracking(true);

        rotate_x += (float)0.02f * dy;

        if(( rotate_x > 360.0f )||( rotate_x < -360.0f ))
        {
            rotate_x = 0.0f;
        }

        rotate_y += (float)0.02f * dx;

        if(( rotate_y > 360.0f )||( rotate_y < -360.0f ))
        {
            rotate_y = 0.0f;
        }

        updateGL();
    }
}
