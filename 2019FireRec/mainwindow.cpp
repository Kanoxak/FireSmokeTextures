#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>
#include <iostream>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core/utility.hpp>

#include "slic.h"

#include <QDebug>

#include "glvisual.h"

/// on_runSeries()
std::ofstream stdof_NFeatures;
std::string stdof_NFields;
std::string tab ("\t");

using namespace cv;
using namespace std;

QString qs_fileName;

std::ofstream stdof_Radius;
std::ofstream stdof_CovRadius;
std::ofstream stdof_Coh;
std::ofstream stdof_CovCoh;

std::ofstream stdof_mopPx;
std::ofstream stdof_mopPy;
std::ofstream stdof_mopQx;
std::ofstream stdof_mopQy;

std::ofstream LUTof_R;
std::ofstream LUTof_G;
std::ofstream LUTof_B;

int iAr_RadMag[21];
int iAr_CohMag[101];

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(&qt_playVDO, SIGNAL(timeout()), this, SLOT(on_play()));

    i_Hmax = 255;
    i_Hmin = 0;
    i_Smax = 255;
    i_Smin = 0;
    i_Vmax = 255;
    i_Vmin = 0;
    i_sb = 13;
    d_vt = 0.1;

    i_vDotR = 255;
    i_vDotG = 0;
    i_vDotB = 0;
    i_vLineR = 0;
    i_vLineG = 255;
    i_vLineB = 0;

    qt_runSeries = new QTimer(this);
    i_row = 0;
    i_col = 0;
    on_drawSeries();
    connect(ui->obckb_runseries, SIGNAL (released()), this, SLOT(runSeries()));
    connect(qt_runSeries, SIGNAL(timeout()), this, SLOT(on_runSeries()));

}

MainWindow::~MainWindow()
{
    delete ui;
}

QString MainWindow::substr_filename(QString qsPathName){
    int iDotPos=0, iSlashPos=0;
    QString qsFileName;
    for (int i = 0; i < qsPathName.size(); ++i) {
        if (qsPathName.at(i) == QChar('/'))
            if (i > iSlashPos)
                iSlashPos = i;

        if (qsPathName.at(i) == QChar('.'))
            if (i > iDotPos)
                iDotPos = i;
    }
    qsFileName = qsPathName.mid(iSlashPos+1, (qsPathName.size()-iSlashPos));

    return qsFileName;
}

QImage MainWindow::on_displayAlpha(cv::Mat mframe)
{
    QImage qframe = QImage(
                        reinterpret_cast<unsigned char*>( mframe.data ),
                        mframe.cols,
                        mframe.rows,
                        QImage::Format_Alpha8
                    );
    return qframe;
}

QImage MainWindow::on_displayRGB(cv::Mat mframe)
{
    QImage qframe = QImage(
                        reinterpret_cast<unsigned char*>( mframe.data ),
                        mframe.cols,
                        mframe.rows,
                        QImage::Format_RGB888
                    ).rgbSwapped();
    return qframe;
}

void MainWindow::on_pushVDOLoad_clicked()
{
    int iFrame = 0;
    double dFramerate = 0.0;
    int iFramewidth = 0;
    int iFrameheight = 0;

    qs_PathName = QFileDialog::getOpenFileName(
                this,
                tr("Open Target Video"),
                "/Users/Phoenix/Documents/QTProjects/videos",
                tr("Video Files (*.avi *.mp4 *.mov *.mpeg *mpg)"),
                nullptr,
                //nullptr//QFileDialog::Options() | QFileDialog::DontUseNativeDialog
                QFileDialog::DontUseSheet |
                QFileDialog::DontUseCustomDirectoryIcons |
                QFileDialog::ReadOnly |
                QFileDialog::HideNameFilterDetails |
                QFileDialog::DontUseNativeDialog);
    cout << "qsPathName: " << qs_PathName.toStdString() << "\n";
    QString qsFileName = qs_fileName = substr_filename( qs_PathName );
    cout << "qs_fileName: " << qsFileName.toStdString() << "\n";

    cv::VideoCapture cap(qs_PathName.toStdString());
    if(!cap.isOpened()){
        cout << "Error opening video stream or file" << endl;
    } else {
        i_frameTotal = int(cap.get(cv::CAP_PROP_FRAME_COUNT));
        i_frame = iFrame = int(cap.get(cv::CAP_PROP_POS_FRAMES));

        dFramerate = cap.get(cv::CAP_PROP_FPS);
        iFramewidth = int(cap.get(3));
        iFrameheight = int(cap.get(4));
        cout << "Successfully opening video stream or file" << endl;
        cout << "Total frame: " << std::to_string(i_frameTotal) << endl;
        cout << "At frame: " << std::to_string(iFrame) << endl;

        cout << "Frame Rate(fps): " << std::to_string(dFramerate) << endl;
        cout << "Frame Resolution: " << std::to_string(iFramewidth) << "x" << std::to_string(iFrameheight) << endl;
    }


    cv::Mat frame;
    /// Capture frame-by-frame
    cap >> frame;

    cv::Mat rgbframe;
    /// rgbSwapped() resulted the image as done by cvtColor() below
    //  cv::cvtColor(frame, rgbframe, CV_BGR2RGB);

    QImage qframe =  on_displayRGB(frame);

    ui->labelVDOPanel->setPixmap(QPixmap::fromImage(qframe));
    ui->labelVDOPanel->setScaledContents(true);
    ui->labelFilename->setText("Filename: " + qs_PathName);
    ui->labelFPS->setText("FPS: " + QString::number(dFramerate));
    ui->labelVDOresolution->setText("Resolution: " + QString::number(iFramewidth) + "x" + QString::number(iFrameheight));
    ui->labelPanelresolution->setText("Panel: " + QString::number(ui->labelVDOPanel->width()) + "x" + QString::number(ui->labelVDOPanel->height()));
    ui->oblbl_frameTotal->setText(QString::number(i_frameTotal));
    ui->oblbl_frame->setText(QString::number(i_frame));

    cap.release();
    i_frame = 0;
    ui->pushPlayVDO->setEnabled(true);
}

QImage Mat2QImage(cv::Mat const& src) {
     cv::Mat temp; // make the same cv::Mat
     cv::cvtColor(src, temp,CV_BGR2RGB); // cvtColor Makes a copt, that what i need
     QImage dest((uchar*) temp.data, temp.cols, temp.rows, temp.step, QImage::Format_RGB888);
     QImage dest2(dest);
     dest2.detach(); // enforce deep copy
     return dest2;
}

QImage MainWindow::on_colorChasing(cv::Mat mframe)
{

    cv::resize(mframe, mframe, cv::Size(ui->labelVDOPanel_9->width(), ui->labelVDOPanel_9->height()));
    cv::cvtColor(mframe, mframe, CV_BGR2HSV);
    cv::Mat mmask = mframe;

    /// Red color cv::Scalar(160, 100, 100) ~ cv::Scalar(179, 255, 255)

    cv::inRange(mframe, cv::Scalar(i_Hmin, i_Smin, i_Vmin), cv::Scalar(i_Hmax, i_Smax, i_Vmax), mmask);
    ui->oblbl_Scalar->setText("cv::Scalar("+ QString::number(i_Hmin) + ","
                                           + QString::number(i_Smin) + ","
                                           + QString::number(i_Vmin) + ")~cv::Scalar("
                                           + QString::number(i_Hmax) + ","
                                           + QString::number(i_Smax) + ","
                                           + QString::number(i_Vmax) + ")");

    QImage qimmasked = on_displayAlpha(mmask);
    ui->labelVDOPanel_9->setPixmap(QPixmap::fromImage(qimmasked));

    cv::cvtColor(mframe, mframe, CV_HSV2RGB);

    cv::Mat mchase = mframe;
    cv::cvtColor(mmask, mmask, CV_GRAY2RGB);

    cv::bitwise_and(mframe, mmask, mchase);

    cv::cvtColor(mchase, mchase, CV_RGB2BGR);
    QImage qimchased = on_displayRGB(mchase);
    return qimchased;
}

double MainWindow::normalizeCoh(double dCoh, double dLbound) {

    // Sept. 6, 2016
    if ((dLbound < dCoh) && (dCoh < 100))
        return 100*(3*pow(((dCoh-dLbound)/(100-dLbound)), 2) - 2*pow(((dCoh-dLbound)/(100-dLbound)), 3));
    else if (dCoh <= dLbound)
        return 0;
    else
        return 100;
}

double absDouble(double d) {
    if (d < 0)
        return -d;
    else
        return d;
}

double minDouble(double a, double b) {
    if (a < b) return a;
    else return b;
}

double MainWindow::coh4Adjacents(double dAngle, const cv::Point2f& pA, const cv::Point2f& pB, const cv::Point2f& pC, double dVt) {

    double dCohA = 0;
    double dCohB = 0;
    double dCohC = 0;

    /// +----+----+
    /// |////| cA |
    /// +----+----+
    /// | cB | cC |
    /// +----+----+

    //int iNeighbors=0;

    if ( (pA.y > 0) || (pA.x > 0) ) {
        double dRadiusPA = sqrt((double) pA.y * (double) pA.y + (double) pA.x * (double) pA.x);
        if (dRadiusPA >= dVt) {
            double dAngleA = atan2((double) pA.y, (double) pA.x);
            if ((dAngleA < 0.0) && (dAngleA > -CV_PI)) {
                dAngleA += (2* CV_PI);
            }
            double dMeanA = ((dAngle + dAngleA)/2.00);
            double dDiff1A = absDouble(dAngle - dMeanA);
            double dDiff2A = absDouble(dAngle - (CV_PI - dMeanA)); // Sept 7, 2016
            //dCohA = absDouble( 100.0-((100.0*minDouble(dDiff1A, dDiff2A))/(CV_PI/2)) );
            dCohA = ( 1-(2.0*minDouble(dDiff1A, dDiff2A)) / CV_PI) *100; // Sept 6, 2016

            //iNeighbors++;

        }
    } else {
        dCohA = 0;
    }

    if ( (pB.y > 0) || (pB.x > 0) ) {
        double dRadiusPB = sqrt((double) pB.y * (double) pB.y + (double) pB.x * (double) pB.x);
        if (dRadiusPB >= dVt) {
            double dAngleB = atan2((double) pB.y, (double) pB.x);
            if ((dAngleB < 0.0) && (dAngleB > -CV_PI)) {
                dAngleB += (2* CV_PI);
            }
            double dMeanB = ((dAngle + dAngleB)/2.00);
            double dDiff1B = absDouble(dAngle - dMeanB);
            double dDiff2B = absDouble(dAngle - (CV_PI - dMeanB)); // Sept 7, 2016
            //dCohB = absDouble( 100.0-((100.0*minDouble(dDiff1B, dDiff2B))/(CV_PI/2)) );
            dCohB = ( 1-(2.0*minDouble(dDiff1B, dDiff2B)) / CV_PI) *100; // Sept 6, 2016

            //iNeighbors++;
        }
    } else {
        dCohB = 0;
    }

    if ( (pC.y > 0) || (pC.x > 0) ) {
        double dRadiusPC = sqrt((double) pC.y * (double) pC.y + (double) pC.x * (double) pC.x);
        if (dRadiusPC >= dVt) {
            double dAngleC = atan2((double) pC.y, (double) pC.x);
            if ((dAngleC < 0.0) && (dAngleC > -CV_PI)) {
                dAngleC += (2* CV_PI);
            }
            double dMeanC = ((dAngle + dAngleC)/2.00);
            double dDiff1C = absDouble(dAngle - dMeanC);
            double dDiff2C = absDouble(dAngle - (CV_PI - dMeanC)); // Sept 7, 2016
            //dCohC = absDouble( 100.0-((100.0*minDouble(dDiff1C, dDiff2C))/(CV_PI/2)) );
            dCohC = ( 1-(2.0*minDouble(dDiff1C, dDiff2C)) / CV_PI) *100; // Sept 6, 2016

            //iNeighbors++;
        }
    } else {
        dCohC = 0;
    }

    ///Sept. 20, 2016 test 2
    /*
    double mean = (dCohA+dCohB+dCohC)/3;
    double CCor = (dCohA+dCohB+dCohC) * mean;
    double SSQ  = (dCohA*dCohA)+(dCohB*dCohB)+(dCohC*dCohC);
    return NormalizeCoh((CCor/SSQ)*mean, 70);
    */

    //return ( (dCohA+dCohB+dCohC)/3 ); // Sept 6, 2016
    double dCohAvg = 0; //June. 6, 2019
    if (( dCohA > 0 ) && ( dCohB > 0 ) && ( dCohC > 0 ))
        dCohAvg = (dCohA + dCohB + dCohC)/3;
    else if (( 0 == dCohA ) && ( dCohB > 0 ) && ( dCohC > 0 ))
        dCohAvg = (dCohB + dCohC)/2;
    else if (( 0 == dCohB ) && ( dCohA > 0 ) && ( dCohC > 0 ))
        dCohAvg = (dCohA + dCohC)/2;
    else if (( 0 == dCohC ) && ( dCohA > 0 ) && ( dCohB > 0 ))
        dCohAvg = (dCohA + dCohB)/2;
    else if (( 0 == dCohA ) && ( 0 == dCohB ) && ( dCohC > 0 ))
        dCohAvg = dCohC;
    else if (( 0 == dCohB ) && ( 0 == dCohC ) && ( dCohA > 0 ))
        dCohAvg = dCohA;
    else if (( 0 == dCohC ) && ( 0 == dCohA ) && ( dCohB > 0 ))
        dCohAvg = dCohB;

    return normalizeCoh( dCohAvg, 70 ); //Sept. 7, 2016


    //return NormalizeCoh( (dCohA+dCohB+dCohC)/iNeighbors, 70); //Sept. 20, 2016 //if iNeighbors <= 3
}

/*
void MainWindow::drawOpticalFlow(const cv::Mat mVector,
                                 cv::Mat mVectorFlow,
                                 int iStep,
                                 const cv::Scalar& scColorLine,
                                 const cv::Scalar& scColorDot,
                                 double dVt)
{
    double dRadius = 0.00;

    QImage img = Mat2QImage(mVectorFlow);

    for(int y = 0; y < mVector.rows; y += iStep)
        for(int x = 0; x < mVector.cols; x += iStep) {
            const cv::Point2f& pFxy = mVector.at<cv::Point2f>(y, x);

            if ((pFxy.x > 0) || (pFxy.y > 0)) {
                int xx = x + pFxy.x;
                int yy = y + pFxy.y;
                dRadius = sqrt(((double) pFxy.y * (double) pFxy.y) + ((double) pFxy.x * (double) pFxy.x));
                if ( dRadius >= dVt ) {
                    circle(mVectorFlow, cvPoint(x, y), 1, scColorDot, 1);
                    line(mVectorFlow, cv::Point(x, y), cv::Point(xx, yy), scColorLine);
                } //if ( dRadius >= dVt ) {
            } // if ((pFxy.x != 0) || (pFxy.y != 0)) {
        } // for(int x = 0; x < mVector.cols; x += iStep) {
    // for(int y = 0; y < mVector.rows; y += iStep) {
}
*/

MainWindow::VectorAttribute MainWindow::drawOpticalFlow(const cv::Mat mVector,
                                 cv::Mat mVectorFlow,
                                 int iStep,
                                 const cv::Scalar& scColorLine,
                                 const cv::Scalar& scColorDot,
                                 double dVt)
{

    VectorAttribute va;
    va.dAngle = 0.00;
    va.sinAngle = 0.00;
    va.cosAngle = 0.00;
    va.dRadius = 0.00;

    double dRadius = 0.00;
    double dDegree = 0.00;
    double dCRadius = 0.00;
    int iCountF = 0;

    double dSin = 0.00, dCos = 0.00;
    double dCoh = 0.00;
    double dSSQ = 0.00;

    for(int jjj=0; jjj<13; jjj++){
        va.daRadius[jjj] = 0;
        if (jjj < 5) {
            va.da5Radius[jjj][0] = 0;
            va.da5Radius[jjj][1] = 0;
        }
    }

    for(int jjj=0; jjj<360; jjj++){
        va.daDegree[jjj] = 0;
        if (jjj < 5) {
            va.da5Degree[jjj][0] = 0;
            va.da5Degree[jjj][1] = 0;
        }
    }

    for(int jjj=0; jjj<256; jjj++){
        va.iaR[jjj] = 0;
        va.iaG[jjj] = 0;
        va.iaB[jjj] = 0;
        if (jjj < 5) {
            va.ia5R[jjj][0] = 0;
            va.ia5G[jjj][0] = 0;
            va.ia5B[jjj][0] = 0;
            va.ia5R[jjj][1] = 0;
            va.ia5G[jjj][1] = 0;
            va.ia5B[jjj][1] = 0;
        }
    }

    ///Feb 16, 2018
    int iRGBtotal=0;
    va.iRMean=0;
    va.iBMean=0;
    va.iGMean=0;

    int iRadius = 0;

    //change to use vector (structure)
    stdof_mopPx.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_pPx.txt");
    stdof_mopPy.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_pPy.txt");
    stdof_mopQx.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_pQx.txt");
    stdof_mopQy.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_pQy.txt");

    if (ui->obckb_SB->isChecked()) {
        LUTof_R.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_R.txt");
        LUTof_G.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_G.txt");
        LUTof_B.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_B.txt");
        //LUTof_Label.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_LUT_Label.txt");
    }

    QImage img = Mat2QImage(mVectorFlow);

    for(int y = 0; y < mVector.rows; y += iStep)
        for(int x = 0; x < mVector.cols; x += iStep) {
            const cv::Point2f& pFxy = mVector.at<cv::Point2f>(y, x);

            if ((pFxy.x > 0) || (pFxy.y > 0)) {
                int xx = x + pFxy.x;
                int yy = y + pFxy.y;
                dRadius = sqrt(((double) pFxy.y * (double) pFxy.y) + ((double) pFxy.x * (double) pFxy.x));
                if ( dRadius >= dVt && ((x < ui->labelVDOPanel_8->width()-5) && (y < ui->labelVDOPanel_8->height()-5) && (xx < ui->labelVDOPanel_8->width()-10) && (yy < ui->labelVDOPanel_8->height()-10)) ) {

                ///if ( dRadius >= dVt ) {
                    circle(mVectorFlow, cvPoint(x, y), 1, scColorDot, 1);
                    line(mVectorFlow, cv::Point(x, y), cv::Point(xx, yy), scColorLine);

                    // Angle in radian not degree
                    dDegree = atan2((double) pFxy.y, (double) pFxy.x);

                    //Sept. 12, 2016.
                    /// Here to store (pFxy, dRadius, pDegree)
                    //stdof_mopix << std::to_string(x) + "|" + std::to_string(y) + "%"
                    //              + std::to_string(dDegree) + "="
                    //              + std::to_string(dRadius) + "\n";

                    stdof_mopPx << std::to_string(x) + "\n";
                    stdof_mopPy << std::to_string(y) + "\n";
                    stdof_mopQx << std::to_string(pFxy.x) + "\n";
                    stdof_mopQy << std::to_string(pFxy.y) + "\n";

                    dSin += sin(dDegree);
                    dCos += cos(dDegree);
                    dCRadius += dRadius;

                    /// +----+----+
                    /// |////| cA |
                    /// +----+----+
                    /// | cB | cC |
                    /// +----+----+

                    const cv::Point2f& pFxyA = mVector.at<cv::Point2f>(y, x+1);
                    const cv::Point2f& pFxyB = mVector.at<cv::Point2f>(y+1, x);
                    const cv::Point2f& pFxyC = mVector.at<cv::Point2f>(y+1, x+1);

                    //Sept. 7, 2016
                    if ((pFxyA.y > 0) || (pFxyA.x > 0) || (pFxyB.y > 0) || (pFxyB.x > 0) || (pFxyC.y > 0) || (pFxyC.x > 0) ) {
                        double dCohIndex = coh4Adjacents( dDegree, pFxyA, pFxyB, pFxyC, dVt );

                        dCoh += dCohIndex;
                        dSSQ += dCohIndex*dCohIndex;
                    }

                    /// Created for testing, no UI for this yet.
                    /*
                    /// +----+----+----+
                    /// | c2 | c3 | c4 |
                    /// +----+----+----+
                    /// | c1 | Ci |////|
                    /// +----+----+----+
                    /// |////|////|////|
                    /// +----+----+----+
                    if ((x > 0) && (y > 0)) {
                        const cv::Point2f& pFxy1 = mFlow.at<cv::Point2f>(y  , x-1);
                        const cv::Point2f& pFxy2 = mFlow.at<cv::Point2f>(y-1, x-1);
                        const cv::Point2f& pFxy3 = mFlow.at<cv::Point2f>(y-1, x  );
                        const cv::Point2f& pFxy4 = mFlow.at<cv::Point2f>(y-1, x+1);

                        dCoh5 += Coh5Adjacents( dDegree, pFxy1, pFxy2, pFxy3, pFxy4, dVt );

                        iCount5++;
                    }
                    */

                    if (ui->obckb_runseries->isChecked()) {
                        for(int q = y; q < y+iStep; q++)
                            for(int p = x; p < x+iStep; p++) {
                              //  QImage img = Mat2QImage(mMappedFlow);
                                int r, g, b;

                                QColor pxl(img.pixel(p, q));

                                //QRgb tempColorRgb = img.pixel(p, q);
                                //QColor tempColor(tempColorRgb);

                                //printf("RGB components of the pixel selected: %d %d %d\n", tempColor.red(), tempColor.green(), tempColor.blue());

                                r = pxl.red();
                                g = pxl.green();
                                b = pxl.blue();

                                /// write R G B pixel by pixel into a .txt file
                                LUTof_R << std::to_string(r) + "\n";
                                LUTof_G << std::to_string(g) + "\n";
                                LUTof_B << std::to_string(b) + "\n";
                                //LUTof_Label << "0\n";

                                ///Feb 16, 2018
                                iRGBtotal++;
                                va.iRMean+=r;
                                va.iGMean+=g;
                                va.iBMean+=b;

                                va.iaR[r]++;
                                va.iaG[g]++;
                                va.iaB[b]++;

                            }
                    } // if (chk_3DVisual->isChecked()) {

                    iCountF++;

                    if ((dRadius >= 0.0) && (dRadius <= 0.5)) {
                        iRadius = 0;
                    } else if ((dRadius > 0.5) && (dRadius <= 1.5)) {
                        iRadius = 1;
                    } else if ((dRadius > 1.5) && (dRadius <= 2.5)) {
                        iRadius = 2;
                    } else if ((dRadius > 2.5) && (dRadius <= 3.5)) {
                        iRadius = 3;
                    } else if ((dRadius > 3.5) && (dRadius <= 4.5)) {
                        iRadius = 4;
                    } else if ((dRadius > 4.5) && (dRadius <= 5.5)) {
                        iRadius = 5;
                    } else if ((dRadius > 5.5) && (dRadius <= 6.5)) {
                        iRadius = 6;
                    } else if ((dRadius > 6.5) && (dRadius <= 7.5)) {
                        iRadius = 7;
                    } else if ((dRadius > 7.5) && (dRadius <= 8.5)) {
                        iRadius = 8;
                    } else if ((dRadius > 8.5) && (dRadius <= 9.5)) {
                        iRadius = 9;
                    } else if ((dRadius > 9.5) && (dRadius <=10.5)) {
                        iRadius =10;
                    } else if ((dRadius >10.5) && (dRadius <=11.5)) {
                        iRadius = 11;
                    } else if (dRadius > 11.5) {
                        iRadius = 12;
                    }

                    va.daRadius[iRadius]++;

                    va.dAngle = atan2(sin(dDegree), cos(dDegree)) * (double)(180.0/CV_PI);
                    if ((va.dAngle < 0.0) && (va.dAngle > -180.0)) {
                        va.dAngle = 360.0 + va.dAngle;
                    }
                    int angle = (int)va.dAngle;
                    if (angle > 359) angle = 0;
                    va.daDegree[angle]++;

                } //if ( dRadius >= dVt ) {
            } // if ((pFxy.x != 0) || (pFxy.y != 0)) {
        } // for(int x = 0; x < mVector.cols; x += iStep) {
    // for(int y = 0; y < mVector.rows; y += iStep) {
    /// Find the most top 5 frequently found R G B Angle Radius
    ///

    double maxim = 0;
    va.ia5R[0][0] = 0;
    va.ia5R[0][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if(va.iaR[iii] > maxim){
            maxim = va.iaR[iii];
            va.ia5R[0][0] = iii;
            va.ia5R[0][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5R[1][0] = 0;
    va.ia5R[1][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaR[iii] >= maxim) && (iii != va.ia5R[0][0])) {
            maxim = va.iaR[iii];
            va.ia5R[1][0] = iii;
            va.ia5R[1][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5R[2][0] = 0;
    va.ia5R[2][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaR[iii] >= maxim) && (iii != va.ia5R[0][0]) && (iii != va.ia5R[1][0])){
            maxim = va.iaR[iii];
            va.ia5R[2][0] = iii;
            va.ia5R[2][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5R[3][0] = 0;
    va.ia5R[3][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaR[iii] >= maxim) && (iii != va.ia5R[0][0]) && (iii != va.ia5R[1][0]) && (iii != va.ia5R[2][0])){
            maxim = va.iaR[iii];
            va.ia5R[3][0] = iii;
            va.ia5R[3][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5R[4][0] = 0;
    va.ia5R[4][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaR[iii] >= maxim) && (iii != va.ia5R[0][0]) && (iii != va.ia5R[1][0]) && (iii != va.ia5R[2][0]) && (iii != va.ia5R[3][0])){
            maxim = va.iaR[iii];
            va.ia5R[4][0] = iii;
            va.ia5R[4][1] = maxim;
        }
    }

    maxim = 0;
    va.ia5G[0][0] = 0;
    va.ia5G[0][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if(va.iaG[iii] > maxim){
            maxim = va.iaG[iii];
            va.ia5G[0][0] = iii;
            va.ia5G[0][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5G[1][0] = 0;
    va.ia5G[1][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaG[iii] >= maxim) && (iii != va.ia5G[0][0])){
            maxim = va.iaG[iii];
            va.ia5G[1][0] = iii;
            va.ia5G[1][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5G[2][0] = 0;
    va.ia5G[2][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaG[iii] >= maxim) && (iii != va.ia5G[0][0]) && (iii != va.ia5G[1][0])){
            maxim = va.iaG[iii];
            va.ia5G[2][0] = iii;
            va.ia5G[2][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5G[3][0] = 0;
    va.ia5G[3][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaG[iii] >= maxim) && (iii != va.ia5G[0][0]) && (iii != va.ia5G[1][0]) && (iii != va.ia5G[2][0])){
            maxim = va.iaG[iii];
            va.ia5G[3][0] = iii;
            va.ia5G[3][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5G[4][0] = 0;
    va.ia5G[4][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaG[iii] >= maxim) && (iii != va.ia5G[0][0]) && (iii != va.ia5G[1][0]) && (iii != va.ia5G[2][0]) && (iii != va.ia5G[3][0])){
            maxim = va.iaG[iii];
            va.ia5G[4][0] = iii;
            va.ia5G[4][1] = maxim;
        }
    }

    maxim = 0;
    va.ia5B[0][0] = 0;
    va.ia5B[0][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if(va.iaB[iii] > maxim){
            maxim = va.iaB[iii];
            va.ia5B[0][0] = iii;
            va.ia5B[0][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5B[1][0] = 0;
    va.ia5B[1][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaB[iii] >= maxim) && (iii != va.ia5B[0][0])){
            maxim = va.iaB[iii];
            va.ia5B[1][0] = iii;
            va.ia5B[1][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5B[2][0] = 0;
    va.ia5B[2][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaB[iii] >= maxim) && (iii != va.ia5B[0][0]) && (iii != va.ia5B[1][0])){
            maxim = va.iaB[iii];
            va.ia5B[2][0] = iii;
            va.ia5B[2][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5B[3][0] = 0;
    va.ia5B[3][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaB[iii] >= maxim) && (iii != va.ia5B[0][0]) && (iii != va.ia5B[1][0]) && (iii != va.ia5B[2][0])){
            maxim = va.iaB[iii];
            va.ia5B[3][0] = iii;
            va.ia5B[3][1] = maxim;
        }
    }
    maxim = 0;
    va.ia5B[4][0] = 0;
    va.ia5B[4][1] = 0;
    for(int iii = 0; iii < 256; iii++){
        if((va.iaB[iii] >= maxim) && (iii != va.ia5B[0][0]) && (iii != va.ia5B[1][0]) && (iii != va.ia5B[2][0]) && (iii != va.ia5B[3][0])){
            maxim = va.iaB[iii];
            va.ia5B[4][0] = iii;
            va.ia5B[4][1] = maxim;
        }
    }

    maxim = 0;
    va.da5Degree[0][0] = 0;
    va.da5Degree[0][1] = 0;
    for(int iii = 0; iii < 360; iii++){
        if(va.daDegree[iii] > maxim){
            maxim = va.daDegree[iii];
            va.da5Degree[0][0] = iii;
            va.da5Degree[0][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Degree[1][0] = 0;
    va.da5Degree[1][1] = 0;
    for(int iii = 0; iii < 360; iii++){
        if((va.daDegree[iii] >= maxim) && (iii != va.da5Degree[0][0])){
            maxim = va.daDegree[iii];
            va.da5Degree[1][0] = iii;
            va.da5Degree[1][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Degree[2][0] = 0;
    va.da5Degree[2][1] = 0;
    for(int iii = 0; iii < 360; iii++){
        if((va.daDegree[iii] >= maxim) && (iii != va.da5Degree[0][0]) && (iii != va.da5Degree[1][0])){
            maxim = va.daDegree[iii];
            va.da5Degree[2][0] = iii;
            va.da5Degree[2][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Degree[3][0] = 0;
    va.da5Degree[3][1] = 0;
    for(int iii = 0; iii < 360; iii++){
        if((va.daDegree[iii] >= maxim) && (iii != va.da5Degree[0][0]) && (iii != va.da5Degree[1][0]) && (iii != va.da5Degree[2][0])){
            maxim = va.daDegree[iii];
            va.da5Degree[3][0] = iii;
            va.da5Degree[3][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Degree[4][0] = 0;
    va.da5Degree[4][1] = 0;
    for(int iii = 0; iii < 360; iii++){
        if((va.daDegree[iii] >= maxim) && (iii != va.da5Degree[0][0]) && (iii != va.da5Degree[1][0]) && (iii != va.da5Degree[2][0]) && (iii != va.da5Degree[3][0])){
            maxim = va.daDegree[iii];
            va.da5Degree[4][0] = iii;
            va.da5Degree[4][1] = maxim;
        }
    }

    maxim = 0;
    va.da5Radius[0][0] = 0;
    va.da5Radius[0][1] = 0;
    for(int iii = 0; iii < 13; iii++){
        if(va.daRadius[iii] > maxim){
            maxim = va.daRadius[iii];
            va.da5Radius[0][0] = iii;
            va.da5Radius[0][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Radius[1][0] = 0;
    va.da5Radius[1][1] = 0;
    for(int iii = 0; iii < 13; iii++){
        if((va.daRadius[iii] > maxim) && (iii != va.da5Radius[0][0])){
            maxim = va.daRadius[iii];
            va.da5Radius[1][0] = iii;
            va.da5Radius[1][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Radius[2][0] = 0;
    va.da5Radius[2][1] = 0;
    for(int iii = 0; iii < 13; iii++){
        if((va.daRadius[iii] > maxim) && (iii != va.da5Radius[0][0]) && (iii != va.da5Radius[1][0])){
            maxim = va.daRadius[iii];
            va.da5Radius[2][0] = iii;
            va.da5Radius[2][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Radius[3][0] = 0;
    va.da5Radius[3][1] = 0;
    for(int iii = 0; iii < 13; iii++){
        if((va.daRadius[iii] > maxim) && (iii != va.da5Radius[0][0]) && (iii != va.da5Radius[1][0]) && (iii != va.da5Radius[2][0])){
            maxim = va.daRadius[iii];
            va.da5Radius[3][0] = iii;
            va.da5Radius[3][1] = maxim;
        }
    }
    maxim = 0;
    va.da5Radius[4][0] = 0;
    va.da5Radius[4][1] = 0;
    for(int iii = 0; iii < 13; iii++){
        if((va.daRadius[iii] > maxim) && (iii != va.da5Radius[0][0]) && (iii != va.da5Radius[1][0]) && (iii != va.da5Radius[2][0]) && (iii != va.da5Radius[3][0])){
            maxim = va.daRadius[iii];
            va.da5Radius[4][0] = iii;
            va.da5Radius[4][1] = maxim;
        }
    }

    if (iCountF > 0) {
        va.dAngle = atan2((dSin/iCountF), (dCos/iCountF)) * (double)(180.0/CV_PI);
        va.sinAngle = dSin/iCountF;
        va.cosAngle = dCos/iCountF;
        va.dRadius = dCRadius/iCountF;
        ///va.dCohIndex = dCoh/iCountF;
        va.dCohIndex = ((dCoh * (dCoh/iCountF)) / dSSQ)*100;
        if (isnan(va.dCohIndex)) va.dCohIndex = 0.000000;
    } else {
        va.dAngle = 0;
        va.sinAngle = 0;
        va.cosAngle = 0;
        va.dRadius = 0;
        va.dCohIndex = 0;
    }
    va.iVecCount = iCountF;

    /// after got dCoh if we want to do more weight using NCCor for motion coherent index on frame
    /// 1. Mean ==> dCoh/iCountF
    /// 2. C = Sum x Mean ==> Sum == dCoh
    /// 3. SSQ = Sum of dCoh^2
    /// 4. NCCor = C/SSQ ==> va.dCohIndex

    va.iRGBTotal = iRGBtotal;
    qDebug() << "va.dCohIndex: " << va.dCohIndex;



    if ((va.dAngle < 0.0) && (va.dAngle > -180.0)) {
        va.dAngle = 360.0 + va.dAngle;
        //va.dAngle = sin(va.dAngle);
    }
    //First passed | first pyramid

    qDebug() << "va.dRadius: " << va.dRadius;
    qDebug() << "dCRadius: " << dCRadius;
    qDebug() << "iCountF: " << iCountF;


    //va.dCohIndex = dCoh5/iCount5;

    if (ui->obckb_runseries->isChecked()) {
        ///Jan 22, 2020
        if(iRGBtotal > 0){
            ///Feb 16, 2018
            va.iRMean/=iRGBtotal;
            va.iGMean/=iRGBtotal;
            va.iBMean/=iRGBtotal;
        }else{
            va.iRMean = 0;
            va.iGMean = 0;
            va.iBMean = 0;
        }
        LUTof_R.close();
        LUTof_G.close();
        LUTof_B.close();
    }

    stdof_mopPx.close();
    stdof_mopPy.close();
    stdof_mopQx.close();
    stdof_mopQy.close();

    return va;
}

QImage MainWindow::drawSuperPixels(cv::Mat mframe) {

    cv::resize(mframe, mframe, cv::Size(ui->labelVDOPanel_11->width(), ui->labelVDOPanel_11->height()));
    cv::Mat mclone = mframe;

    /* Load the image and convert to Lab colour space. */
    cv::cvtColor(mframe, mclone, CV_BGR2Lab);

     /* Yield the number of superpixels and weight-factors from the user. */
    int w = mframe.cols, h = mframe.rows;
    int nr_superpixels = 100;

    /* The colour (nc) parameters. */
    int nc = 64; // 256 / 8

    double step = sqrt(( w * h) / (double) nr_superpixels);

    /* Perform the SLIC superpixel algorithm. */
    Slic slic;
    /* The step size per cluster, and the colour (nc) and distance (ns)
     * parameters. */
    slic.generate_superpixels(mclone, step, nc);
    slic.create_connectivity(mclone);

    cv::cvtColor(mframe, mframe, CV_Lab2BGR);
    /* Display the contours and show the result. */
    slic.display_contours(mframe, CV_RGB(0, 255, 0));
    slic.display_center_grid(mframe, CV_RGB(255, 255, 0));

    QImage qimSuperPix = on_displayRGB(mframe);
    return qimSuperPix;
}

cv::Mat3b QImage2Mat(const QImage &src) {
    unsigned int height = src.height();
    unsigned int width = src.width();

    cv::Mat3b dest(height, width);
    for (unsigned int y = 0; y < height; ++y) {
        cv::Vec3b *destrow = dest[y];
        for (unsigned int x = 0; x < width; ++x) {
            QRgb pxl = src.pixel(x, y);
            destrow[x] = cv::Vec3b(qBlue(pxl), qGreen(pxl), qRed(pxl));
        }
    }
    return dest;
}

void MainWindow::on_play() {
    cv::VideoCapture cap(qs_PathName.toStdString());
    cv::VideoCapture prevcap(qs_PathName.toStdString());

    if (!cap.isOpened()) {
        cout << "Error opening video stream or file" << endl;
    } else {
        if (0 == i_frame) {
            i_frameTotal = int(cap.get(cv::CAP_PROP_FRAME_COUNT));
            double iFramerate = cap.get(cv::CAP_PROP_FPS);
            int iFramewidth = int(cap.get(3));
            int iFrameheight = int(cap.get(4));
            cout << "Successfully opening video stream or file" << endl;
            cout << "Total frame: " << std::to_string(i_frameTotal) << endl;
            cout << "At frame: " << std::to_string(i_frame) << endl;
            cout << "Frame Rate(fps): " << std::to_string(iFramerate) << endl;
            cout << "Frame Resolution: " << std::to_string(iFramewidth) << "x" << std::to_string(iFrameheight) << endl;
            cout << "=========================================================" << endl;
        }
    }

    /// Capture frame-by-frame
    cap.set(cv::CAP_PROP_POS_FRAMES, i_frame);
    cap.set(cv::CAP_PROP_BUFFERSIZE, 5);
    cap >> m_frame;

    //  cout << "At frame: " << std::to_string(i_frame) << endl;

    if (i_frame > 0) {
        prevcap.set(cv::CAP_PROP_POS_FRAMES, i_frame-1);
        prevcap.set(cv::CAP_PROP_BUFFERSIZE, 5);
        prevcap >> m_prevframe;

        cv::Mat mdiff, mframe, mprevframe;

        cv::cvtColor(m_frame, mframe, CV_BGR2GRAY);
        cv::cvtColor(m_prevframe, mprevframe, CV_BGR2GRAY);
        QImage qframe_2 = on_displayAlpha(mprevframe);
        ui->labelVDOPanel_2->clear();
        ui->labelVDOPanel_2->setPixmap(QPixmap::fromImage(qframe_2));
        ui->labelVDOPanel_2->setScaledContents(true);
        ui->labelVDOPanel_2->repaint();

        if(ui->obckb_absdiff->isChecked()){
            cv::absdiff(mframe, mprevframe, mdiff);
            QImage qframe_3 = on_displayAlpha(mdiff);
            ui->labelVDOPanel_3->clear();
            ui->labelVDOPanel_3->setPixmap(QPixmap::fromImage(qframe_3));
            ui->labelVDOPanel_3->setScaledContents(true);
            ui->labelVDOPanel_3->repaint();
        }

        if(ui->obckb_bitwiseXOR->isChecked()){
            cv::bitwise_xor(mframe, mprevframe, mdiff);
            QImage qframe_4 = on_displayAlpha(mdiff);
            ui->labelVDOPanel_4->clear();
            ui->labelVDOPanel_4->setPixmap(QPixmap::fromImage(qframe_4));
            ui->labelVDOPanel_4->setScaledContents(true);
            ui->labelVDOPanel_4->repaint();
        }

        if(ui->obckb_subtract->isChecked()){
            cv::subtract(mframe, mprevframe, mdiff, cv::Mat(), mdiff.type());
            QImage qframe_5 = on_displayAlpha(mdiff);
            ui->labelVDOPanel_5->clear();
            ui->labelVDOPanel_5->setPixmap(QPixmap::fromImage(qframe_5));
            ui->labelVDOPanel_5->setScaledContents(true);
            ui->labelVDOPanel_5->repaint();
        }
        /// November 12, 2019 Notes:
        /// This threshold below used the difference of a pair of frames calculated by subtract() method.
        if(ui->obckb_threshold->isChecked()){
            if(ui->rbt_absdiff->isChecked())
                cv::absdiff(mframe, mprevframe, mdiff);
            else if(ui->rbt_bitwiseXOR->isChecked())
                cv::bitwise_xor(mframe, mprevframe, mdiff);
            else
                cv::subtract(mframe, mprevframe, mdiff, cv::Mat(), mdiff.type());

            cv::threshold(mdiff, mdiff, 10, 255, cv::THRESH_BINARY);
            QImage qframe_6 = on_displayAlpha(mdiff);
            ui->labelVDOPanel_6->clear();
            ui->labelVDOPanel_6->setPixmap(QPixmap::fromImage(qframe_6));
            ui->labelVDOPanel_6->setScaledContents(true);
            ui->labelVDOPanel_6->repaint();
        }
        /// November 12, 2019 Notes:
        /// For the FC... videos, the continity may drop at pair of frames when index as 1, 7, 13, ... increasing with +6.
        /// Those key-frames reported that the number of non-zero pixels were less than around <100.
        /// After carefully looking at the cctv... video, they both reported that the non-zero pixels were sometimes
        /// less than hundred. It depends on more or less motion in front of the camera.
        /// So, we need to consider this challenge and find a solution for one compatible coding.
        /// For example, a solution we can easily utilize is to set the threshold of displaying by < 70.
        /// It means if less motion the system will not detect.
        /*
        vector<cv::Point> whitee_pixels;
        cv::findNonZero(mdiff, whitee_pixels);
        cout << i_frame << " " << whitee_pixels.size() << endl;
        */

        /// November 14, 2019 Notes:
        /// At this moment, if the new video loading happened it will be error at calc optflowgd.cpp:1114.
        /// The question is how m_framesequence has the attribute "data", what does it mean?
        /// Run one error one, new loading is interupted the process of calc optical flow "Farneback".
        ///
        /// Example:
        /// cv::Mat img1 = cv::imread(argv[1], cv::IMREAD_GRAYSCALE);
        /// cv::Mat img2 = cv::imread(argv[2], cv::IMREAD_GRAYSCALE);
        /// cv::Mat res;
        ///
        /// cv::calcOpticalFlowFarneback(img1, img2, res, .4, 1, 12, 2, 8, 1.2, 0);
        ///
        /// cv:imshow("cat", res);
        /// cv::waitKey(0);
        ///
        /// Add the following line before cv::calcOpticalFlowFarneback().
        /// cv::resize(img2, img2, img1.size());
        ///               

        cv::Mat mframeTOI = m_prevframe;
        cv::resize(mframeTOI, mframeTOI, cv::Size(ui->labelVDOPanel_8->width(), ui->labelVDOPanel_8->height()));
        cv::Mat mframeTOIGray = m_prevframe;
        cv::resize(mframeTOIGray, mframeTOIGray, cv::Size(ui->labelVDOPanel_8->width(), ui->labelVDOPanel_8->height()));        

        QImage qimmask = on_colorChasing(m_prevframe);

        /// cv::Mat mframeHist = m_prevframe;
        /// cv::resize(mframeHist, mframeHist, cv::Size(ui->labelVDOPanel_10->width(), ui->labelVDOPanel_10->height()));
        /// cv::cvtColor(mframeHist, mframeHist, CV_RGB2GRAY);
        /// cv::equalizeHist( mframeHist, mframeHist );
        /// cv::cvtColor(mframeHist, mframeHist, CV_GRAY2RGB);
        /// if converted GRAY --> RGB don't forget to change on_displayRGB( )
        /// QImage qimHist = on_displayAlpha(mframeHist);
        /// ui->labelVDOPanel_10->setPixmap(QPixmap::fromImage(qimHist));
        /// ui->labelVDOPanel_10->setScaledContents(true);

        if (ui->obckb_SPix->isChecked()) {
            QImage qimSuperPixels = drawSuperPixels(QImage2Mat(qimmask));
            ui->labelVDOPanel_11->setPixmap(QPixmap::fromImage(qimSuperPixels));
            ui->labelVDOPanel_11->setScaledContents(true);
        }


        mframeTOI = QImage2Mat(qimmask);

        /// If using mask just change from mframeTOI --> mmask.
        cv::cvtColor(mframeTOI, mframeTOIGray, CV_RGB2GRAY);

        /// For testing using the equalized hist grayscale frame
        /// cv::cvtColor(mframeHist, mframeHist, CV_RGB2GRAY);
        /// mframeTOIGray = mframeHist;

        if ( m_framesequence.data ) {
            /// In the paper needs to describe this part clearly and more.

            /// Comment of the different of his and my parameters
            /*
            void cvCalcOpticalFlowFarneback( const CvArr* prev, const CvArr* next, CvArr* flow,
                                             double     dpyr_scale,
                                             int        ilevels,
                                             int        iwinsize,
                                             int        iiterations,
                                             int        ipoly_n,
                                             double     dpoly_sigma,
                                             int        iflags)
            cv::calcOpticalFlowFarneback( m_framesequence, mframeTOIGray, m_vector,
                                          0.4,
                                          1,
                                          12,
                                          32,
                                          8,
                                          1.2,
                                          0);
            cv::calcOpticalFlowFarneback( m_framesequence, mframeTOIGray, m_vector,
                                          0.5,
                                          3,
                                          15,
                                          3,
                                          5,
                                          1.2,
                                          0);
            */
            cv::calcOpticalFlowFarneback( m_framesequence, mframeTOIGray, m_vector,
                                          0.5,
                                          3,
                                          15,
                                          3,
                                          5,
                                          1.2,
                                          0);
            va_Res = MainWindow::drawOpticalFlow(
                                    m_vector,
                                    mframeTOI,
                                    i_sb,
                                    CV_RGB(i_vLineR, i_vLineG, i_vLineB),
                                    CV_RGB(i_vDotR, i_vDotG, i_vDotB),
                                    d_vt
                                );

            QImage qframe = on_displayRGB(mframeTOI);
            ui->labelVDOPanel_8->clear();
            ui->labelVDOPanel_8->setPixmap(QPixmap::fromImage(qframe));
            ui->labelVDOPanel_8->setScaledContents(true);
            ui->labelVDOPanel_8->repaint();

            if (ui->obckb_runseries->isChecked() && va_Res.dRadius > d_vt) {

                if(i_begin < i_end) {
                    d_AvgRadius += dAr_Rad[i_begin] = va_Res.dRadius;
                    ///d_AvgAngle += dAr_Ang[i_begin] = va_Res.dAngle;
                    d_AvgCoh += dAr_Coh[i_begin] = va_Res.dCohIndex;
                    ///d_AvgMCoh += dAr_MCoh[i_begin] = d_MeanCoh;

                    i_begin++;
                } else if(i_begin == i_end && !b_out) { //one time do

                    b_out = true;
                    d_AvgRadius /= i_end;
                    ///d_AvgAngle /= i_end;
                    d_AvgCoh /= i_end;
                    ///d_AvgMCoh /= i_end;

                    for (int i = 0; i < i_end; i++) {

                        d_StdRad += (dAr_Rad[i] - d_AvgRadius) * (dAr_Rad[i] - d_AvgRadius);
                        ///d_StdAng += (dAr_Ang[i] - d_AvgAngle) * (dAr_Ang[i] - d_AvgAngle);
                        d_StdCoh += (dAr_Coh[i] - d_AvgCoh) * (dAr_Coh[i] - d_AvgCoh);
                        ///d_StdMCoh += (dAr_MCoh[i] - d_AvgMCoh) * (dAr_MCoh[i] - d_AvgMCoh);

                    }

                    d_StdRad = sqrt(double(1.00/(double)i_end) * d_StdRad);
                    ///d_StdAng = sqrt(double(1.00/(double)i_end) * d_StdAng);
                    d_StdCoh = sqrt(double(1.00/(double)i_end) * d_StdCoh);
                    ///d_StdMCoh = sqrt(double(1.00/(double)i_end) * d_StdMCoh);

                }

            }

        } std::swap( m_framesequence, mframeTOIGray );
    }

    QImage qframe = QImage(reinterpret_cast<const unsigned char*>(m_frame.data),
                    m_frame.cols,
                    m_frame.rows,
                    QImage::Format_RGB888).rgbSwapped();

    ui->labelVDOPanel->clear();
    ui->labelVDOPanel->setPixmap(QPixmap::fromImage(qframe));
    ui->labelVDOPanel->setScaledContents(true);
    ui->labelVDOPanel->repaint();

    i_frame++;

    int i_frameSld = 100*i_frame/i_frameTotal;
    ui->obhsl_frame->setValue(i_frameSld);
    ui->oblbl_frameTotal->setText(QString::number(i_frameTotal));
    ui->oblbl_frame->setText(QString::number(i_frame));

    if ( int(cap.get(cv::CAP_PROP_FRAME_COUNT)) == i_frame) {
        qt_playVDO.stop();
        cap.release();
    }
}

void MainWindow::on_pushPlayVDO_clicked()
{
    qt_playVDO.start();
    i_frame = 0;
    on_play();
}

void MainWindow::on_obhsl_frame_sliderMoved(int position)
{
    // 760 bar length
    i_frame = cvCeil((position * i_frameTotal)/100);
    ui->oblbl_frame->setText(QString::number(i_frame));
}

void MainWindow::on_obhsl_frame_sliderPressed()
{
    i_frame = cvCeil(( ui->obhsl_frame->value() * i_frameTotal ) /100);
    ui->oblbl_frame->setText(QString::number(i_frame));
}

void MainWindow::on_obhsl_frame_sliderReleased()
{
    i_frame = cvCeil(( ui->obhsl_frame->value() * i_frameTotal ) /100);
    ui->oblbl_frame->setText(QString::number(i_frame));
}

void MainWindow::on_obvsl_Hmin_sliderMoved(int position)
{
    i_Hmin = position;
    ui->oblbl_Hmin->setText(QString::number(i_Hmin));
}

void MainWindow::on_obvsl_Hmin_sliderPressed()
{
    i_Hmin = ui->obvsl_Hmin->value();
    ui->oblbl_Hmin->setText(QString::number(i_Hmin));
}

void MainWindow::on_obvsl_Hmin_sliderReleased()
{
    i_Hmin = ui->obvsl_Hmin->value();
    ui->oblbl_Hmin->setText(QString::number(i_Hmin));
}

void MainWindow::on_obvsl_Smin_sliderMoved(int position)
{
    i_Smin = position;
    ui->oblbl_Smin->setText(QString::number(i_Smin));
}

void MainWindow::on_obvsl_Smin_sliderPressed()
{
    i_Smin = ui->obvsl_Smin->value();
    ui->oblbl_Smin->setText(QString::number(i_Smin));
}

void MainWindow::on_obvsl_Smin_sliderReleased()
{
    i_Smin = ui->obvsl_Smin->value();
    ui->oblbl_Smin->setText(QString::number(i_Smin));
}

void MainWindow::on_obvsl_Vmin_sliderMoved(int position)
{
    i_Vmin = position;
    ui->oblbl_Vmin->setText(QString::number(i_Vmin));
}

void MainWindow::on_obvsl_Vmin_sliderPressed()
{
    i_Vmin = ui->obvsl_Vmin->value();
    ui->oblbl_Vmin->setText(QString::number(i_Vmin));
}

void MainWindow::on_obvsl_Vmin_sliderReleased()
{
    i_Vmin = ui->obvsl_Vmin->value();
    ui->oblbl_Vmin->setText(QString::number(i_Vmin));
}

void MainWindow::on_obvsl_Hmax_sliderMoved(int position)
{
    i_Hmax = position;
    ui->oblbl_Hmax->setText(QString::number(i_Hmax));
}

void MainWindow::on_obvsl_Hmax_sliderPressed()
{
    i_Hmax = ui->obvsl_Hmax->value();
    ui->oblbl_Hmax->setText(QString::number(i_Hmax));
}

void MainWindow::on_obvsl_Hmax_sliderReleased()
{
    i_Hmax = ui->obvsl_Hmax->value();
    ui->oblbl_Hmax->setText(QString::number(i_Hmax));
}

void MainWindow::on_obvsl_Smax_sliderMoved(int position)
{
    i_Smax = position;
    ui->oblbl_Smax->setText(QString::number(i_Smax));
}

void MainWindow::on_obvsl_Smax_sliderPressed()
{
    i_Smax = ui->obvsl_Smax->value();
    ui->oblbl_Smax->setText(QString::number(i_Smax));
}

void MainWindow::on_obvsl_Smax_sliderReleased()
{
    i_Smax = ui->obvsl_Smax->value();
    ui->oblbl_Smax->setText(QString::number(i_Smax));
}

void MainWindow::on_obvsl_Vmax_sliderMoved(int position)
{
    i_Vmax = position;
    ui->oblbl_Vmax->setText(QString::number(i_Vmax));
}

void MainWindow::on_obvsl_Vmax_sliderPressed()
{
    i_Vmax = ui->obvsl_Vmax->value();
    ui->oblbl_Vmax->setText(QString::number(i_Vmax));
}

void MainWindow::on_obvsl_Vmax_sliderReleased()
{
    i_Vmax = ui->obvsl_Vmax->value();
    ui->oblbl_Vmax->setText(QString::number(i_Vmax));
}

void MainWindow::on_obhsl_sb_sliderMoved(int position)
{
    /// i_sb length is between 3-21 step 2
    int isb = ((position/2)*2)+1;
    i_sb = isb;
    ui->oblbl_sb->setText(QString::number(i_sb));
}

void MainWindow::on_obhsl_sb_sliderPressed()
{
     int isb = ((ui->obhsl_sb->value()/2)*2)+1;
     i_sb = isb;
     ui->oblbl_sb->setText(QString::number(i_sb));
}

void MainWindow::on_obhsl_sb_sliderReleased()
{
    int isb = ((ui->obhsl_sb->value()/2)*2)+1;
    i_sb = isb;
    ui->oblbl_sb->setText(QString::number(i_sb));
}

void MainWindow::on_obhsl_vt_sliderMoved(int position)
{
    /// d_vt length is between 0-20 step 1 ~ 0.0-2.0
    //  double dvt = ( (double)position/10.0 );
    d_vt = ( (double)position/10.0 );
    ui->oblbl_vt->setText(QString::number(d_vt));
}

void MainWindow::on_obhsl_vt_sliderPressed()
{
    d_vt = ( (double)ui->obhsl_vt->value()/10.0 );
    ui->oblbl_vt->setText(QString::number(d_vt));
}

void MainWindow::on_obhsl_vt_sliderReleased()
{
    d_vt = ( (double)ui->obhsl_vt->value()/10.0 );
    ui->oblbl_vt->setText(QString::number(d_vt));
}

void MainWindow::on_drawSeries()
{
    CvFont cvfFont;
    double dHScale = 0.3;
    double dVScale = 0.3;
    int iLineWidth = 1;
    cvInitFont(&cvfFont, CV_FONT_ITALIC|CV_FONT_NORMAL, dHScale, dVScale, 0, iLineWidth);

    int gray = 220;

    /// Radius Magnitude
    ///cv::Mat mSeriesRM = Mat::zeros(ui->oblbl_line01RM->height(), ui->oblbl_line01RM->width(), CV_8UC3);
    ///mSeriesRM.setTo(Scalar::all(255));
    cv::Mat mSeriesRM = Mat(ui->oblbl_line01RM->height(), ui->oblbl_line01RM->width(), CV_8UC3, cvScalar(255, 255, 255));
    for (int y=13; y<228; y+=20) {
        line(mSeriesRM, cv::Point(30, y), cv::Point(433, y), CV_RGB(gray, gray, gray), 1);
    }
    for (int x=33; x<440; x+=20) {
        line(mSeriesRM,  cv::Point(x, 13), cv::Point(x, 217), CV_RGB(gray, gray, gray), 1);
    }
    ///QImage qiMonitorUL = on_displayRGB(mSeriesUL); // because it is not image but it is drawing
    QImage qiMonitorRM = Mat2QImage(mSeriesRM);
    ui->oblbl_line01RM->setPixmap(QPixmap::fromImage(qiMonitorRM));

    /// Radius
    cv::Mat mSeriesR = Mat(ui->oblbl_line01R->height(), ui->oblbl_line01R->width(), CV_8UC3, cvScalar(255, 255, 255));
    for (int y=13; y<128; y+=20) {
        cv::line(mSeriesR,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
    }
    for (int x=33; x<400; x+=20) {
        cv::line(mSeriesR,  cv::Point(x, 13), cv::Point(x, 117), CV_RGB(gray, gray, gray), 1);
    }
    ///-- Radius Covariance Stationary --///
    for (int y=138; y<228; y+=20) {
        cv::line(mSeriesR,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
        if (178==y) {
            cv::line(mSeriesR,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
        }
    }
    for (int x=33; x<400; x+=20) {
        cv::line(mSeriesR,  cv::Point(x, 134), cv::Point(x, 218), CV_RGB(gray, gray, gray), 1);
    }
    QImage qiMonitorR = Mat2QImage(mSeriesR);
    ui->oblbl_line01R->setPixmap(QPixmap::fromImage(qiMonitorR));

    /// Coherent Magnitude 2019
    cv::Mat mSeriesCM = Mat(ui->oblbl_line01CM->height(), ui->oblbl_line01CM->width(), CV_8UC3, cvScalar(255, 255, 255));
    for (int y=13; y<228; y+=20) {
        cv::line(mSeriesCM,  cv::Point(30, y), cv::Point(433, y), CV_RGB(gray, gray, gray), 1);
    }
    for (int x=33; x<440; x+=20) {
        cv::line(mSeriesCM,  cv::Point(x, 13), cv::Point(x, 217), CV_RGB(gray, gray, gray), 1);
    }
    QImage qiMonitorCM = Mat2QImage(mSeriesCM);
    ui->oblbl_line01CM->setPixmap(QPixmap::fromImage(qiMonitorCM));

    /// Coherent index
    cv::Mat mSeriesC = Mat(ui->oblbl_line01C->height(), ui->oblbl_line01C->width(), CV_8UC3, cvScalar(255, 255, 255));
    for (int y=13; y<128; y+=20) {
        cv::line(mSeriesC,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
    }
    for (int x=33; x<400; x+=20) {
        cv::line(mSeriesC,  cv::Point(x, 13), cv::Point(x, 117), CV_RGB(gray, gray, gray), 1);
    }
    ///-- Coherence Index Covariance Stationary --///
    for (int y=138; y<228; y+=20) {
        cv::line(mSeriesC,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
        if (178==y) {
            cv::line(mSeriesC,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
        }
    }
    for (int x=33; x<400; x+=20) {
        cv::line(mSeriesC,  cv::Point(x, 134), cv::Point(x, 218), CV_RGB(gray, gray, gray), 1);
    }
    QImage qiMonitorC = Mat2QImage(mSeriesC);
    ui->oblbl_line01C->setPixmap(QPixmap::fromImage(qiMonitorC));

    /// Sin Cos
    cv::Mat mSeriesCS = Mat(ui->oblbl_line01CS->height(), ui->oblbl_line01CS->width(), CV_8UC3, cvScalar(255, 255, 255));
    for (int y=13; y<228; y+=20) {

        if (113!=y) {
            cv::line(mSeriesCS,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
        }
        if (53==y) {
            cv::line(mSeriesCS,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
        }
        if (173==y) {
            cv::line(mSeriesCS,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
        }
    }
    for (int x=33; x<400; x+=20) {
        cv::line(mSeriesCS,  cv::Point(x, 13), cv::Point(x, 107), CV_RGB(gray, gray, gray), 1);
    }
    for (int x=33; x<400; x+=20) {
        cv::line(mSeriesCS,  cv::Point(x, 120), cv::Point(x, 213), CV_RGB(gray, gray, gray), 1);
    }
    QImage qiMonitorCS = Mat2QImage(mSeriesCS);
    ui->oblbl_line01CS->setPixmap(QPixmap::fromImage(qiMonitorCS));
}


int iCount;
void MainWindow::runSeries()
{
        if (ui->obckb_runseries->isChecked()) {
            ///on_drawSeries();
            b_out = false;
            i_begin = 0;
            i_end = 30;
            d_AvgRadius = 0;
            ///d_AvgAngle = 0;
            d_AvgCoh = 0;
            ///d_AvgMCoh = 0;
            d_StdRad = 0;
            ///d_StdAng = 0;
            d_StdCoh = 0;
            ///d_StdMCoh = 0;


            b_firstReleased = false;
            on_runSeries();
            qt_runSeries->start();
            iCount = 0;

        } else {
            b_out = true;
            i_begin = i_end = 0;
            on_runSeries();
            b_firstReleased = false;
            qt_runSeries->stop();
        }
}

cv::Point cvpRadMag;
cv::Mat mSeriesRM;
int iFirstRadMag;

cv::Point cvpRadius;
cv::Mat mSeriesR;
int iFirstRadius;

cv::Point cvpCohIndex;
cv::Mat mSeriesC;
int iFirstCohIndex;

cv::Point cvpCohMag;
cv::Mat mSeriesCM;
int iFirstCohMag;

cv::Point cvpCos;
cv::Point cvpSin;
cv::Mat mSeriesCS;
int iFirstCos;
int iFirstSin;

///cv::Point cvpVarAngle;
///int iFirstVarAngle;
cv::Point cvpVarRadius;
int iFirstVarRadius;
cv::Point cvpVarCohIndex;
int iFirstVarCohIndex;
///cv::Point cvpVarMeanCoh;
///int iFirstVarMeanCoh;

double dtmp_Rad = 0;
double dtmp_Coh = 0;

///int iAr_RadMag[21];
///int iAr_CohMag[101];

void MainWindow::on_runSeries()
{
    // Dec. 29, 2017
    double d_varRad=0;
    double d_varCoh=0;

    // Jan. 25, 2018
    ///double d_varAng=0;
    double d_pearsoncorcoef;

    qDebug() << "iCount++ :: " << iCount;
    //cout << "@qs_fileName: " << qs_fileName.toStdString() << "\n";


    if(0 == iCount) {
        for (int iAr=0; iAr<=20; iAr++) {
            iAr_RadMag[iAr] = 0;
            iAr_CohMag[iAr] = 0;
        }

        for (int iAc=21; iAc<=100; iAc++){
            iAr_CohMag[iAc] = 0;
        }
    }

    /// if 30 first frames have passed
    /// one time do
    if ( i_begin == i_end ) {
        if (!b_firstReleased) {

            mSeriesR = Mat(ui->oblbl_line01R->height(), ui->oblbl_line01R->width(), CV_8UC3, cvScalar(255, 255, 255));
            mSeriesRM = Mat(ui->oblbl_line01RM->height(), ui->oblbl_line01RM->width(), CV_8UC3, cvScalar(255, 255, 255));
            mSeriesC = Mat(ui->oblbl_line01C->height(), ui->oblbl_line01C->width(), CV_8UC3, cvScalar(255, 255, 255));
            mSeriesCM = Mat(ui->oblbl_line01CM->height(), ui->oblbl_line01CM->width(), CV_8UC3, cvScalar(255, 255, 255));
            mSeriesCS = Mat(ui->oblbl_line01CS->height(), ui->oblbl_line01CS->width(), CV_8UC3, cvScalar(255, 255, 255));

            //cout << "@va_Res.dRadius: " << va_Res.dRadius << "\n";

            b_firstReleased = true;
            cvpCohMag.x = cvpCohIndex.x = cvpRadius.x = cvpRadMag.x = iCount = 0;
            cvpRadMag.y = iFirstRadMag = (int) round(va_Res.dRadius * 200 / 360);
            cvpRadius.y = iFirstRadius = (int) round(va_Res.dRadius * 100 / 20);
            cvpCohIndex.y = iFirstCohIndex = (int) round(va_Res.dCohIndex);
            cvpCohMag.y = iFirstCohMag = (int) round(va_Res.dCohIndex * 200 / 360);
            cvpCos.y = iFirstCos = (int) round(va_Res.cosAngle * 10 * 4.0);
            cvpSin.y = iFirstSin = (int) round(va_Res.sinAngle * 10 * 4.0);


            cvpVarCohIndex.x = cvpVarRadius.x = 0; ///cvpVarMeanCoh.x =  ... cvpVarAngle.x =
            ///cvpVarAngle.y = iFirstVarAngle = (int) round(((va_Res.dAngle - d_AvgAngle) / d_StdAng) * 4.0);
            if (0.0 == d_StdRad){
                cvpVarRadius.y = iFirstVarRadius = 0;
            }else{
                cvpVarRadius.y = iFirstVarRadius = (int) round(((va_Res.dRadius - d_AvgRadius) / d_StdRad) * 4.0);
            }

            if (0.0 == d_StdCoh){
                cvpVarCohIndex.y = iFirstVarCohIndex = 0;
            }else{
                cvpVarCohIndex.y = iFirstVarCohIndex = (int) round(((va_Res.dCohIndex - d_AvgCoh) / d_StdCoh) * 4.0);
            }
            ///cvpVarMeanCoh.y = iFirstVarMeanCoh = (int) round(((d_MeanCoh - d_AvgMCoh) / d_StdMCoh) * 4.0);


            CvFont cvfFont;
            double dHScale = 0.3;
            double dVScale = 0.3;
            int iLineWidth = 1;
            cvInitFont(&cvfFont, CV_FONT_ITALIC|CV_FONT_NORMAL, dHScale, dVScale, 0, iLineWidth);
            int gray = 240;

            /// Radius Magnitude
            ///cv::Mat mSeriesRM = Mat::zeros(ui->oblbl_line01RM->height(), ui->oblbl_line01RM->width(), CV_8UC3);
            ///mSeriesRM.setTo(Scalar::all(255));
            //cv::Mat mSeriesRM = Mat(ui->oblbl_line01RM->height(), ui->oblbl_line01RM->width(), CV_8UC3, cvScalar(255, 255, 255));
            for (int y=13; y<228; y+=20) {
                line(mSeriesRM, cv::Point(30, y), cv::Point(433, y), CV_RGB(gray, gray, gray), 1);
            }
            for (int x=33; x<440; x+=20) {
                line(mSeriesRM,  cv::Point(x, 13), cv::Point(x, 217), CV_RGB(gray, gray, gray), 1);
            }

            /// Radius
            //cv::Mat mSeriesR = Mat(ui->oblbl_line01R->height(), ui->oblbl_line01R->width(), CV_8UC3, cvScalar(255, 255, 255));
            for (int y=13; y<128; y+=20) {
                cv::line(mSeriesR,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
            }
            for (int x=33; x<400; x+=20) {
                cv::line(mSeriesR,  cv::Point(x, 13), cv::Point(x, 117), CV_RGB(gray, gray, gray), 1);
            }
            ///-- Radius Covariance Stationary --///
            for (int y=138; y<228; y+=20) {
                cv::line(mSeriesR,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
                if (178==y) {
                    cv::line(mSeriesR,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
                }
            }
            for (int x=33; x<400; x+=20) {
                cv::line(mSeriesR,  cv::Point(x, 134), cv::Point(x, 218), CV_RGB(gray, gray, gray), 1);
            }

            /// Coherent Magnitude 2019
            //cv::Mat mSeriesCM = Mat(ui->oblbl_line01CM->height(), ui->oblbl_line01CM->width(), CV_8UC3, cvScalar(255, 255, 255));
            for (int y=13; y<228; y+=20) {
                cv::line(mSeriesCM,  cv::Point(30, y), cv::Point(433, y), CV_RGB(gray, gray, gray), 1);
            }
            for (int x=33; x<440; x+=20) {
                cv::line(mSeriesCM,  cv::Point(x, 13), cv::Point(x, 217), CV_RGB(gray, gray, gray), 1);
            }

            /// Coherent index
            //cv::Mat mSeriesC = Mat(ui->oblbl_line01C->height(), ui->oblbl_line01C->width(), CV_8UC3, cvScalar(255, 255, 255));
            for (int y=13; y<128; y+=20) {
                cv::line(mSeriesC,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
            }
            for (int x=33; x<400; x+=20) {
                cv::line(mSeriesC,  cv::Point(x, 13), cv::Point(x, 117), CV_RGB(gray, gray, gray), 1);
            }

            ///-- Coherence Index Covariance Stationary --///
            for (int y=138; y<228; y+=20) {
                cv::line(mSeriesC,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
                if (178==y) {
                    cv::line(mSeriesC,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
                }
            }
            for (int x=33; x<400; x+=20) {
                cv::line(mSeriesC,  cv::Point(x, 134), cv::Point(x, 218), CV_RGB(gray, gray, gray), 1);
            }

            /// Sin Cos
            //cv::Mat mSeriesCS = Mat(ui->oblbl_line01CS->height(), ui->oblbl_line01CS->width(), CV_8UC3, cvScalar(255, 255, 255));
            for (int y=13; y<228; y+=20) {

                if (113!=y) {
                    cv::line(mSeriesCS,  cv::Point(30, y), cv::Point(393, y), CV_RGB(gray, gray, gray), 1);
                }
                if (53==y) {
                    cv::line(mSeriesCS,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
                }
                if (173==y) {
                    cv::line(mSeriesCS,  cv::Point(30, y), cv::Point(393, y), CV_RGB(205, 170, 0), 1);
                }
            }
            for (int x=33; x<400; x+=20) {
                cv::line(mSeriesCS,  cv::Point(x, 13), cv::Point(x, 107), CV_RGB(gray, gray, gray), 1);
            }
            for (int x=33; x<400; x+=20) {
                cv::line(mSeriesCS,  cv::Point(x, 120), cv::Point(x, 213), CV_RGB(gray, gray, gray), 1);
            }
        } else {
            if (va_Res.dRadius > d_vt ) {

                ///  33
                ///13 +----+----+  20 =  13
                ///   | 20 |    |
                ///   +----+----+  16 =  33
                ///   |    |    |     :
                ///   +----+----+   0 = 113
                ///   rows = 100 pixels
                ///   cols = 300 pixels

                /// ::Radius (0-20)
                cvpRadius.y = (int) round(va_Res.dRadius * 100 / 20);
                qDebug() << "cvpRadius.y: " << cvpRadius.y;
                cv::line(mSeriesR, cv::Point(33+cvpRadius.x, 113-iFirstRadius), cv::Point(33+cvpRadius.x, 113-cvpRadius.y), CV_RGB(12, 16, 213), 1);

                /// ::Coherent Index (0-100)
                cvpCohIndex.y = (int) round(va_Res.dCohIndex);
                cv::line(mSeriesC, cv::Point(33+cvpCohIndex.x, 113-iFirstCohIndex), cv::Point(33+cvpCohIndex.x, 113-cvpCohIndex.y), CV_RGB(13, 17, 172), 1);


                if(iCount <= 360) {
                    int ir = (int) round(va_Res.dRadius);
                    iAr_RadMag[ir] = (int) iAr_RadMag[ir] + 1;
                    /// qDebug() << "iAr_RadMag[" << ir << "]" << iAr_RadMag[ir] << "::" << (int) round(va_Res.dRadius);

                    int ic = (int) round(va_Res.dCohIndex);
                    iAr_CohMag[ic] = (int) iAr_CohMag[ic] + 1;
                    /// qDebug() << "iAr_CohMag[" << ic << "]" << iAr_CohMag[ic] << "::" << (int) round(va_Res.dCohIndex);

                    if(360 == iCount){
                        for (int iAr=0; iAr<=20; iAr++) {
                            cvpRadMag.y = (int) round(iAr_RadMag[iAr] * 200 / 360);
                            if (iAr>0) {
                                cv::line(mSeriesRM, cv::Point(33+cvpRadMag.x, 213-iFirstRadMag), cv::Point(33+cvpRadMag.x+20, 213-cvpRadMag.y), CV_RGB(0, 0, 255), 1);
                                cvpRadMag.x = iAr*20;
                            }
                            iFirstRadMag = cvpRadMag.y;

                            cvpCohMag.y = (int) round(iAr_CohMag[iAr] * 200 / 360);
                            if(iAr>0){
                                cv::line(mSeriesCM, cv::Point(33+cvpCohMag.x, 213-iFirstCohMag), cv::Point(33+cvpCohMag.x+4, 213-cvpCohMag.y), CV_RGB(0, 0, 255), 1);
                                cvpCohMag.x = iAr*4;
                            }
                            iFirstCohMag = cvpCohMag.y;
                        }

                        for(int iAc=21; iAc<=100; iAc++){
                            cvpCohMag.y = (int) round(iAr_CohMag[iAc] * 200 / 360);
                            if(iAc>0){
                                cv::line(mSeriesCM, cv::Point(33+cvpCohMag.x, 213-iFirstCohMag), cv::Point(33+cvpCohMag.x+4, 213-cvpCohMag.y), CV_RGB(0, 0, 255), 1);
                                cvpCohMag.x = iAc*4;
                            }
                            iFirstCohMag = cvpCohMag.y;
                        }
                    }
                    qDebug() << "iAr_RadMag[" << ir << "]" << iAr_RadMag[ir] << "::" << (int) round(va_Res.dRadius);
                    qDebug() << "iAr_CohMag[" << ic << "]" << iAr_CohMag[ic] << "::" << (int) round(va_Res.dCohIndex);
                }

                int y1 = 0;
                int y2 = 0;
                /// Monitor &::COS SIN (-1.0 - 1.0)
                cvpSin.y = (int) round(va_Res.sinAngle * 10.0 * 4.0);

                if(iFirstSin > 0) { // <=0; May. 15, 2019
                    y1 = 53-abs(iFirstSin); // May 15, 2019;
                } else {
                    y1 = 53+abs(iFirstSin);
                }
                if(cvpSin.y > 0) {
                    y2 = 53-abs(cvpSin.y);
                } else {
                    y2 = 53+abs(cvpSin.y);
                }
                if(0 == iCount) y1 = 53;
                if(iCount > 0){
                    cv::line(mSeriesCS, cv::Point(33+cvpSin.x, y1), cv::Point(33+cvpSin.x, y2), CV_RGB(172, 17, 13), 1);
                }

                cvpCos.y = (int) round(va_Res.cosAngle * 10.0 * 4.0);

                if(iFirstCos > 0) { // <=0; May. 15, 2019
                    y1 = 173-abs(iFirstCos); // May 15, 2019;
                } else {
                    y1 = 173+abs(iFirstCos);
                }
                if(cvpCos.y > 0) {
                    y2 = 173-abs(cvpCos.y);
                } else {
                    y2 = 173+abs(cvpCos.y);
                }
                if(0 == iCount) y1 = 173;
                if(iCount > 0){
                    cv::line(mSeriesCS, cv::Point(33+cvpCos.x, y1), cv::Point(33+cvpCos.x, y2), CV_RGB(172, 70, 13), 1);
                }

                cvpCohIndex.x = cvpRadius.x = cvpSin.x = cvpCos.x = iCount;
                iFirstRadius = cvpRadius.y;
                iFirstCohIndex = cvpCohIndex.y;
                iFirstSin = cvpSin.y;
                iFirstCos = cvpCos.y;

                // Dec. 29, 2017
                // Jan. 24, 2020
                if (d_StdRad == 0.000000){
                    d_varRad = 0.000000;
                }else{
                    d_varRad = (double)(( va_Res.dRadius - d_AvgRadius ) / d_StdRad);
                }
                if (d_StdCoh == 0.000000){
                    d_varCoh = 0.000000;
                }else{
                    d_varCoh = (double)(( va_Res.dCohIndex - d_AvgCoh ) / d_StdCoh);
                }
                if (isnan(d_varCoh)) {
                    d_varCoh = 0.000000;
                }


                // Jan. 25, 2018
                ///d_varAng = ( va_Res.dAngle - d_AvgAngle ) / d_StdAng;

                if (( d_StdRad * d_StdCoh) == 0.000000){
                    d_pearsoncorcoef = 0.000000;
                }else{
                    d_pearsoncorcoef = (double)(( d_varRad * d_varCoh ) / ( d_StdRad * d_StdCoh));
                }
                qDebug() << "( d_StdRad * d_StdCoh ) :: ( " << d_StdRad << " * " << d_StdCoh << " )";
                cvpVarRadius.y = (int) round( d_varRad * 4.0 );
                if (cvpVarRadius.y > 47)
                    cvpVarRadius.y = 47;
                if (cvpVarRadius.y <-47)
                    cvpVarRadius.y =-47;

                if(iFirstVarRadius > 0) { // <=0; Sept. 9, 2016
                    y1 = 178-abs(iFirstVarRadius); // Dec 9, 2017; 168 --> 166
                } else {
                    y1 = 178+abs(iFirstVarRadius);
                }
                if(cvpVarRadius.y > 0) { // <=0; Sept. 9, 2016
                    y2 = 178-abs(cvpVarRadius.y);
                } else {
                    y2 = 178+abs(cvpVarRadius.y);
                }
                if(1 == iCount) y1 = 178;
                cv::line(mSeriesR, cv::Point(33+cvpVarRadius.x, y1), cv::Point(33+cvpVarRadius.x, y2), CV_RGB(213, 16, 12), 1);

                cvpVarCohIndex.y = (int) round( d_varCoh * 4.0 );
                if (cvpVarCohIndex.y > 47)
                    cvpVarCohIndex.y = 47;
                if (cvpVarCohIndex.y <-47)
                    cvpVarCohIndex.y =-47;

                if(iFirstVarCohIndex > 0) { // <=0; Sept. 9, 2016
                    y1 = 178-abs(iFirstVarCohIndex); // Dec 9, 2017; 168 --> 166
                } else {
                    y1 = 178+abs(iFirstVarCohIndex);
                }
                if(cvpVarCohIndex.y > 0) { // <=0; Sept. 9, 2016
                    y2 = 178-abs(cvpVarCohIndex.y);
                } else {
                    y2 = 178+abs(cvpVarCohIndex.y);
                }
                if(1 == iCount) y1 = 178;
                cv::line(mSeriesC, cv::Point(33+cvpVarCohIndex.x, y1), cv::Point(33+cvpVarCohIndex.x, y2), CV_RGB(172, 17, 13), 1);

                // Dec 7, 2017
                cvpVarCohIndex.x = cvpVarRadius.x = iCount; ///cvpVarMeanCoh.x = ... cvpVarAngle.x =
                ///iFirstVarAngle = cvpVarAngle.y;
                iFirstVarRadius = cvpVarRadius.y;
                iFirstVarCohIndex = cvpVarCohIndex.y;
                ///iFirstVarMeanCoh = cvpVarMeanCoh.y;
            }
        }

        if((( iCount >= 1 && iCount <= 1000 ) && va_Res.dRadius > d_vt) && (va_Res.iVecCount > 10)){

            stdof_NFeatures.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "NF.csv", std::ios::out | std::ios::app);

            if ( 1 == iCount ) {
                stdof_NFeatures <<  "label\tiVecCoun" + tab +

                                    "dRadius" + tab +
                                    "dCohIndex" + tab +
                                    ///"d_MeanCoh" + tab +

                                    "vardRadius" + tab +
                                    "vardCohIndex" + tab +
                                    ///"vard_MeanCoh" + tab +

                                    "d_varRad" + tab +
                                    "d_varCoh" + tab +
                                    ///"d_varMCoh" + tab +

                                    "da5Radius_v0" + tab +
                                    "da5Radius_v1" + tab +
                                    "da5Radius_v2" + tab +
                                    "da5Radius_v3" + tab +
                                    "da5Radius_v4" + tab +
                                    "da5Radius_f0" + tab +
                                    "da5Radius_f1" + tab +
                                    "da5Radius_f2" + tab +
                                    "da5Radius_f3" + tab +
                                    "da5Radius_f4" + tab +
                                    "da5Radius_%0" + tab +
                                    "da5Radius_%1" + tab +
                                    "da5Radius_%2" + tab +
                                    "da5Radius_%3" + tab +
                                    "da5Radius_%4" + tab +

                                    "dCos" + tab +
                                    "dSin" + tab +

                                    "da5Degree_v0" + tab +
                                    "da5Degree_v1" + tab +
                                    "da5Degree_v2" + tab +
                                    "da5Degree_v3" + tab +
                                    "da5Degree_v4" + tab +
                                    "da5Degree_f0" + tab +
                                    "da5Degree_f1" + tab +
                                    "da5Degree_f2" + tab +
                                    "da5Degree_f3" + tab +
                                    "da5Degree_f4" + tab +
                                    "da5Degree_%0" + tab +
                                    "da5Degree_%1" + tab +
                                    "da5Degree_%2" + tab +
                                    "da5Degree_%3" + tab +
                                    "da5Degree_%4" + tab +

                                    "iRGBTotal" + tab +
                                    "iRMean" + tab +
                                    "iGMean" + tab +
                                    "iBMean" + tab +
                                    "d_luminance" + tab +

                                    "ia5R_v0" + tab +
                                    "ia5R_v1" + tab +
                                    "ia5R_v2" + tab +
                                    "ia5R_v3" + tab +
                                    "ia5R_v4" + tab +
                                    "ia5R_f0" + tab +
                                    "ia5R_f1" + tab +
                                    "ia5R_f2" + tab +
                                    "ia5R_f3" + tab +
                                    "ia5R_f4" + tab +
                                    "ia5R_%0" + tab +
                                    "ia5R_%1" + tab +
                                    "ia5R_%2" + tab +
                                    "ia5R_%3" + tab +
                                    "ia5R_%4" + tab +

                                    "ia5G_v0" + tab +
                                    "ia5G_v1" + tab +
                                    "ia5G_v2" + tab +
                                    "ia5G_v3" + tab +
                                    "ia5G_v4" + tab +
                                    "ia5G_f0" + tab +
                                    "ia5G_f1" + tab +
                                    "ia5G_f2" + tab +
                                    "ia5G_f3" + tab +
                                    "ia5G_f4" + tab +
                                    "ia5G_%0" + tab +
                                    "ia5G_%1" + tab +
                                    "ia5G_%2" + tab +
                                    "ia5G_%3" + tab +
                                    "ia5G_%4" + tab +

                                    "ia5B_v0" + tab +
                                    "ia5B_v1" + tab +
                                    "ia5B_v2" + tab +
                                    "ia5B_v3" + tab +
                                    "ia5B_v4" + tab +
                                    "ia5B_f0" + tab +
                                    "ia5B_f1" + tab +
                                    "ia5B_f2" + tab +
                                    "ia5B_f3" + tab +
                                    "ia5B_f4" + tab +
                                    "ia5B_%0" + tab +
                                    "ia5B_%1" + tab +
                                    "ia5B_%2" + tab +
                                    "ia5B_%3" + tab +
                                    "ia5B_%4" + "\n";

                dtmp_Rad = 0;
                dtmp_Coh = 0;
                ///dtmp_MCoh = 0;

            } // if (iCount == 1)

            stdof_Radius.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_y_Radius.txt", std::ios::out | std::ios::app);
            stdof_CovRadius.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_y_CovRadius.txt", std::ios::out | std::ios::app);
            stdof_Coh.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_y_Coh.txt", std::ios::out | std::ios::app);
            stdof_CovCoh.open("/Users/Phoenix/Documents/QTProjects/outtext/"+ qs_fileName.toStdString() + "_y_CovCoh.txt", std::ios::out | std::ios::app);

            /// In every 30 frames, we will analyse their SVD
            /// using
            /// cv::SVD::compute(YY, S, U, V) or
            /// cv::SVD homographySVD(homography, cv::SVD::FULL_UV); // constructor or
            /// homographySVD(newHomography, cv::SVD::FULL_UV); // operator ()
            /// homographySVD.w.at<double>(0, 0); // access the first singular value
            /// alternatives:
            /// cv::SVD::compute(homography, w); // compute just the singular values
            /// cv::eigen(homography, w);

            /// How to do:
            /// create a matrix 10x30 of those four features:
            /// 1. radius
            /// 2. motion coherence index
            /// 3. covariance of 30 prior frames (should be any 30 consecutive frames) of radius
            /// 4. covariance of 30 prior frames (should be any 30 consecutive frames) of motion coherence index

            // Dec 18, 2017
            /// write y values into a .txt file
            /// 1 is fire
            /// 2 is smoke
            /// 3 is waterfall

            //stdof_4Features << "1 1:" + std::to_string( va_Res.dRadius ) +
            //                   " 2:" + std::to_string( va_Res.dCohIndex ) +
            //                   " 3:" + std::to_string( absDouble( d_varRad ) ) +
            //                   " 4:" + std::to_string( absDouble( d_varCoh ) ) +
            //                   "\n";

            /// Feb 25,2018
            //stdof_6Features << "1 1:" + std::to_string( va_Res.dRadius ) +
            //                   " 2:" + std::to_string( va_Res.dCohIndex ) +
            //                   " 3:" + std::to_string( d_MeanCoh ) +
            //                   " 4:" + std::to_string( absDouble( d_varRad ) ) +
            //                   " 5:" + std::to_string( absDouble( d_varCoh ) ) +
            //                   " 6:" + std::to_string( absDouble( d_varMCoh ) ) +
            //                   "\n";

            /// Feb 17, 2018
            //stdof_7Features << "1 1:" + std::to_string( va_Res.dRadius ) +
            //                   " 2:" + std::to_string( va_Res.dCohIndex ) +
            //                   " 3:" + std::to_string( absDouble( d_varRad ) ) +
            //                   " 4:" + std::to_string( absDouble( d_varCoh ) ) +
            //                   " 5:" + std::to_string( va_Res.iRMean ) +
            //                   " 6:" + std::to_string( va_Res.iGMean ) +
            //                   " 7:" + std::to_string( va_Res.iBMean ) +
            //                   "\n";

            /// The intensity of light emitted from a surface or called luminance
            /// is weighted sum of three corresponding linear intensity values (R G B) as in equation below
            /// I = 0.299R + 0.587G + 0.114B
            /// 0 - 255 channels
            ///
            ///
            /// d_varAng = ( va_Res.dAngle - d_AvgAngle ) / d_StdAng;
            /// d_varMCoh = ( d_MeanCoh - d_AvgMCoh ) / d_StdMCoh;
            ///
            /// if (std::isnan(va_Res.dAngle))

            double d_luminance = 0.299*va_Res.iRMean + 0.587*va_Res.iGMean + 0.114*va_Res.iBMean;

            /// label 0 -- first stage of fire where smoke--> infi
            /// label 1 -- second stage of fire where smoke > fire and smoke--> infi
            /// label 2 -- third stage of fire where smoke < fire and fire--> infi
            /// label 3 -- forth stage of fire where fire--> infi

            /// the first line here contains one label and one attribute
            ///
            double da5Radius[5], da5Degree[5];
            int ia5R[5], ia5G[5], ia5B[5];


            for(int i = 0; i < 5; i++){
                da5Radius[i] = 0;
                da5Degree[i] = 0;
                ia5R[i] = 0;
                ia5G[i] = 0;
                ia5B[i] = 0;
            }

            qDebug() << "(double)va_Res.iVecCount :: " << (double)va_Res.iVecCount;
            if(va_Res.iVecCount > 0){
                da5Radius[0] = ( ((double)(va_Res.da5Radius[0][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Radius[1] = ( ((double)(va_Res.da5Radius[1][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Radius[2] = ( ((double)(va_Res.da5Radius[2][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Radius[3] = ( ((double)(va_Res.da5Radius[3][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Radius[4] = ( ((double)(va_Res.da5Radius[4][1] * 100.00) / (double)va_Res.iVecCount) );

                da5Degree[0] = ( ((double)(va_Res.da5Degree[0][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Degree[1] = ( ((double)(va_Res.da5Degree[1][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Degree[2] = ( ((double)(va_Res.da5Degree[2][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Degree[3] = ( ((double)(va_Res.da5Degree[3][1] * 100.00) / (double)va_Res.iVecCount) );
                da5Degree[4] = ( ((double)(va_Res.da5Degree[4][1] * 100.00) / (double)va_Res.iVecCount) );
            }else{
                da5Radius[0] = 0;
                da5Radius[1] = 0;
                da5Radius[2] = 0;
                da5Radius[3] = 0;
                da5Radius[4] = 0;

                da5Degree[0] = 0;
                da5Degree[1] = 0;
                da5Degree[2] = 0;
                da5Degree[3] = 0;
                da5Degree[4] = 0;
            }
            qDebug() << "(double)va_Res.iRGBTotal :: " << (double)va_Res.iRGBTotal;
            if(va_Res.iRGBTotal > 0){
                ia5R[0] = ( ((double)(va_Res.ia5R[0][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5R[1] = ( ((double)(va_Res.ia5R[1][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5R[2] = ( ((double)(va_Res.ia5R[2][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5R[3] = ( ((double)(va_Res.ia5R[3][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5R[4] = ( ((double)(va_Res.ia5R[4][1] * 100.00) / (double)va_Res.iRGBTotal) );

                ia5G[0] = ( ((double)(va_Res.ia5G[0][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5G[1] = ( ((double)(va_Res.ia5G[1][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5G[2] = ( ((double)(va_Res.ia5G[2][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5G[3] = ( ((double)(va_Res.ia5G[3][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5G[4] = ( ((double)(va_Res.ia5G[4][1] * 100.00) / (double)va_Res.iRGBTotal) );

                ia5B[0] = ( ((double)(va_Res.ia5B[0][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5B[1] = ( ((double)(va_Res.ia5B[1][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5B[2] = ( ((double)(va_Res.ia5B[2][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5B[3] = ( ((double)(va_Res.ia5B[3][1] * 100.00) / (double)va_Res.iRGBTotal) );
                ia5B[4] = ( ((double)(va_Res.ia5B[4][1] * 100.00) / (double)va_Res.iRGBTotal) );
            }else{
                ia5R[0] = 0;
                ia5R[1] = 0;
                ia5R[2] = 0;
                ia5R[3] = 0;
                ia5R[4] = 0;

                ia5G[0] = 0;
                ia5G[1] = 0;
                ia5G[2] = 0;
                ia5G[3] = 0;
                ia5G[4] = 0;

                ia5B[0] = 0;
                ia5B[1] = 0;
                ia5B[2] = 0;
                ia5B[3] = 0;
                ia5B[4] = 0;
            }


            stdof_NFeatures << "1\t" + std::to_string( va_Res.iVecCount ) +

                                "\t" + std::to_string( va_Res.dRadius ) +
                                "\t" + std::to_string( va_Res.dCohIndex ) +
                                ///"\t" + std::to_string( d_MeanCoh ) +

                                "\t" + std::to_string( absDouble(va_Res.dRadius-dtmp_Rad) ) +
                                "\t" + std::to_string( absDouble(va_Res.dCohIndex-dtmp_Coh) ) +
                                ///"\t" + std::to_string( absDouble(d_MeanCoh-dtmp_MCoh) ) +

                                "\t" + std::to_string( d_varRad ) + //absDouble() May 28, 2019
                                "\t" + std::to_string( d_varCoh ) + //absDouble() May 28, 2019
                                ///"\t" + std::to_string( d_varMCoh ) + //absDouble() May 28, 2019

                                "\t" + std::to_string( va_Res.da5Radius[0][0] ) +
                                "\t" + std::to_string( va_Res.da5Radius[1][0] ) +
                                "\t" + std::to_string( va_Res.da5Radius[2][0] ) +
                                "\t" + std::to_string( va_Res.da5Radius[3][0] ) +
                                "\t" + std::to_string( va_Res.da5Radius[4][0] ) +
                                "\t" + std::to_string( va_Res.da5Radius[0][1] ) +
                                "\t" + std::to_string( va_Res.da5Radius[1][1] ) +
                                "\t" + std::to_string( va_Res.da5Radius[2][1] ) +
                                "\t" + std::to_string( va_Res.da5Radius[3][1] ) +
                                "\t" + std::to_string( va_Res.da5Radius[4][1] ) +

                                "\t" + std::to_string( da5Radius[0] ) +
                                "\t" + std::to_string( da5Radius[1] ) +
                                "\t" + std::to_string( da5Radius[2] ) +
                                "\t" + std::to_string( da5Radius[3] ) +
                                "\t" + std::to_string( da5Radius[4] ) +

                                "\t" + std::to_string( va_Res.cosAngle ) +
                                "\t" + std::to_string( va_Res.sinAngle ) +

                                "\t" + std::to_string( va_Res.da5Degree[0][0] ) +
                                "\t" + std::to_string( va_Res.da5Degree[1][0] ) +
                                "\t" + std::to_string( va_Res.da5Degree[2][0] ) +
                                "\t" + std::to_string( va_Res.da5Degree[3][0] ) +
                                "\t" + std::to_string( va_Res.da5Degree[4][0] ) +
                                "\t" + std::to_string( va_Res.da5Degree[0][1] ) +
                                "\t" + std::to_string( va_Res.da5Degree[1][1] ) +
                                "\t" + std::to_string( va_Res.da5Degree[2][1] ) +
                                "\t" + std::to_string( va_Res.da5Degree[3][1] ) +
                                "\t" + std::to_string( va_Res.da5Degree[4][1] ) +
                                "\t" + std::to_string( da5Degree[0] ) +
                                "\t" + std::to_string( da5Degree[1] ) +
                                "\t" + std::to_string( da5Degree[2] ) +
                                "\t" + std::to_string( da5Degree[3] ) +
                                "\t" + std::to_string( da5Degree[4] ) +

                                "\t" + std::to_string( va_Res.iRGBTotal ) +
                                "\t" + std::to_string( va_Res.iRMean ) +
                                "\t" + std::to_string( va_Res.iGMean ) +
                                "\t" + std::to_string( va_Res.iBMean ) +
                                "\t" + std::to_string( d_luminance ) +

                                "\t" + std::to_string( va_Res.ia5R[0][0] ) +
                                "\t" + std::to_string( va_Res.ia5R[1][0] ) +
                                "\t" + std::to_string( va_Res.ia5R[2][0] ) +
                                "\t" + std::to_string( va_Res.ia5R[3][0] ) +
                                "\t" + std::to_string( va_Res.ia5R[4][0] ) +
                                "\t" + std::to_string( va_Res.ia5R[0][1] ) +
                                "\t" + std::to_string( va_Res.ia5R[1][1] ) +
                                "\t" + std::to_string( va_Res.ia5R[2][1] ) +
                                "\t" + std::to_string( va_Res.ia5R[3][1] ) +
                                "\t" + std::to_string( va_Res.ia5R[4][1] ) +
                                "\t" + std::to_string( ia5R[0] ) +
                                "\t" + std::to_string( ia5R[1] ) +
                                "\t" + std::to_string( ia5R[2] ) +
                                "\t" + std::to_string( ia5R[3] ) +
                                "\t" + std::to_string( ia5R[4] ) +

                                "\t" + std::to_string( va_Res.ia5G[0][0] ) +
                                "\t" + std::to_string( va_Res.ia5G[1][0] ) +
                                "\t" + std::to_string( va_Res.ia5G[2][0] ) +
                                "\t" + std::to_string( va_Res.ia5G[3][0] ) +
                                "\t" + std::to_string( va_Res.ia5G[4][0] ) +
                                "\t" + std::to_string( va_Res.ia5G[0][1] ) +
                                "\t" + std::to_string( va_Res.ia5G[1][1] ) +
                                "\t" + std::to_string( va_Res.ia5G[2][1] ) +
                                "\t" + std::to_string( va_Res.ia5G[3][1] ) +
                                "\t" + std::to_string( va_Res.ia5G[4][1] ) +
                                "\t" + std::to_string( ia5G[0] ) +
                                "\t" + std::to_string( ia5G[1] ) +
                                "\t" + std::to_string( ia5G[2] ) +
                                "\t" + std::to_string( ia5G[3] ) +
                                "\t" + std::to_string( ia5G[4] ) +

                                "\t" + std::to_string( va_Res.ia5B[0][0] ) +
                                "\t" + std::to_string( va_Res.ia5B[1][0] ) +
                                "\t" + std::to_string( va_Res.ia5B[2][0] ) +
                                "\t" + std::to_string( va_Res.ia5B[3][0] ) +
                                "\t" + std::to_string( va_Res.ia5B[4][0] ) +
                                "\t" + std::to_string( va_Res.ia5B[0][1] ) +
                                "\t" + std::to_string( va_Res.ia5B[1][1] ) +
                                "\t" + std::to_string( va_Res.ia5B[2][1] ) +
                                "\t" + std::to_string( va_Res.ia5B[3][1] ) +
                                "\t" + std::to_string( va_Res.ia5B[4][1] ) +
                                "\t" + std::to_string( ia5B[0] ) +
                                "\t" + std::to_string( ia5B[1] ) +
                                "\t" + std::to_string( ia5B[2] ) +
                                "\t" + std::to_string( ia5B[3] ) +
                                "\t" + std::to_string( ia5B[4] ) +

                               "\n";

            dtmp_Rad = va_Res.dRadius;
            dtmp_Coh = va_Res.dCohIndex;

            stdof_Radius << std::to_string( va_Res.dRadius ) + "\n";
            stdof_CovRadius << std::to_string( d_varRad ) + "\n";
            stdof_Coh << std::to_string( va_Res.dCohIndex ) + "\n";
            stdof_CovCoh << std::to_string( d_varCoh ) + "\n";

            stdof_NFeatures.close();
            stdof_Radius.close();
            stdof_Coh.close();
            stdof_CovRadius.close();
            stdof_CovCoh.close();

            qDebug() << "stdof_CovCoh :: " << d_varCoh;

        }// if( iCount >= 1 && iCount <= 1000 )

        if( iCount <= 360 ) {
            QImage qiSeriesRM = Mat2QImage(mSeriesRM);
            ui->oblbl_line01RM->setPixmap(QPixmap::fromImage(qiSeriesRM));

            QImage qiSeriesR = Mat2QImage(mSeriesR);
            ui->oblbl_line01R->setPixmap(QPixmap::fromImage(qiSeriesR));

            QImage qiSeriesC = Mat2QImage(mSeriesC);
            ui->oblbl_line01C->setPixmap(QPixmap::fromImage(qiSeriesC));

            QImage qiSeriesCM = Mat2QImage(mSeriesCM);
            ui->oblbl_line01CM->setPixmap(QPixmap::fromImage(qiSeriesCM));

            QImage qiSeriesCS = Mat2QImage(mSeriesCS);
            ui->oblbl_line01CS->setPixmap(QPixmap::fromImage(qiSeriesCS));
        }

        if((va_Res.dRadius > d_vt) && (va_Res.iVecCount > 10)){
            qDebug() << "iCount++ :: " << iCount;
            iCount++;
        }

    }// if (i_Begin == i_End)
}

void MainWindow::on_obckb_runseries_stateChanged(int arg1)
{
    if (ui->obckb_runseries->isChecked() and arg1 == 1) {
        b_out = false;
        i_begin = 0;
        i_end = 30;

        d_AvgRadius = 0;
        ///d_AvgAngle = 0;
        d_AvgCoh = 0;
        ///d_AvgMCoh = 0;
        d_StdRad = 0;
        ///d_StdAng = 0;
        d_StdCoh = 0;
        ///d_StdMCoh = 0;
        on_runSeries();
        qt_runSeries->start();
        iCount = 0;

    } else {
        b_out = true;
        i_begin = i_end = 0;

        on_runSeries();
        b_firstReleased = false;

        qt_runSeries->stop();
    }
}
