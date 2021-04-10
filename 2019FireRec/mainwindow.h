#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtCore/qglobal.h>
#if QT_VERSION >= 0x050000
#include <QtWidgets/QMainWindow>
#else
#include <QtGui/QMainWindow>
#endif

#include <opencv2/core/base.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/core/types.hpp>

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/objdetect.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/videoio.hpp>

#include <opencv2/highgui/highgui_c.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/videoio/videoio_c.h>
#include <QTimer>


/// on_runSeries()
extern int iCount;

extern cv::Point cvpRadMag;
extern cv::Mat mSeriesRM;
extern int iFirstRadMag;

extern cv::Point cvpRadius;
extern cv::Mat mSeriesR;
extern int iFirstRadius;

extern cv::Point cvpCohIndex;
extern cv::Mat mSeriesC;
extern int iFirstCohIndex;

extern cv::Point cvpCohMag;
extern cv::Mat mSeriesCM;
extern int iFirstCohMag;

extern cv::Point cvpVarAngle;
extern int iFirstVarAngle;
extern cv::Point cvpVarRadius;
extern int iFirstVarRadius;
extern cv::Point cvpVarCohIndex;
extern int iFirstVarCohIndex;
extern cv::Point cvpVarMeanCoh;
extern int iFirstVarMeanCoh;

extern cv::Point cvpCos;
extern cv::Point cvpSin;
extern cv::Mat mSeriesCS;
extern int iFirstCos;
extern int iFirstSin;

extern int iCoh;

extern double dtmp_Rad;
extern double dtmp_Coh;

extern QString qs_fileName;

extern std::ofstream stdof_NFeatures;
extern std::string stdof_NFields;
extern std::string tab;

extern std::ofstream stdof_Radius;
extern std::ofstream stdof_CovRadius;
extern std::ofstream stdof_Coh;
extern std::ofstream stdof_CovCoh;

extern std::ofstream stdof_mopPx;
extern std::ofstream stdof_mopPy;
extern std::ofstream stdof_mopQx;
extern std::ofstream stdof_mopQy;

extern std::ofstream LUTof_R;
extern std::ofstream LUTof_G;
extern std::ofstream LUTof_B;

extern int iAr_RadMag[21];
extern int iAr_CohMag[101];


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow* ui;
    struct VectorAttribute {
        double dAngle;
        double sinAngle;
        double cosAngle;
        double dRadius;
        double dCohIndex;
        int iVecCount;
        int iRMean;
        int iGMean;
        int iBMean;
        int iaR[256];
        int iaG[256];
        int iaB[256];
        double daDegree[360];

        double daRadius[13];
        int ia5R[5][2];
        int ia5G[5][2];
        int ia5B[5][2];
        double da5Degree[5][2];
        double da5Radius[5][2];
        int iRGBTotal;
    };
    VectorAttribute va_Res;
    double d_MeanCoh;

private slots:
    QString substr_filename(QString);
    void    on_pushVDOLoad_clicked();
    void    on_pushPlayVDO_clicked();
    void    on_play();
    QImage  on_displayRGB(cv::Mat);
    QImage  on_displayAlpha(cv::Mat);

    VectorAttribute    drawOpticalFlow(const cv::Mat, cv::Mat, int, const cv::Scalar&, const cv::Scalar&, double);
    QImage  on_colorChasing(cv::Mat);
    QImage  drawSuperPixels(cv::Mat);

    void    on_obhsl_frame_sliderMoved(int position);
    void    on_obvsl_Hmin_sliderMoved(int position);
    void    on_obvsl_Smin_sliderMoved(int position);
    void    on_obvsl_Vmin_sliderMoved(int position);
    void    on_obvsl_Hmax_sliderMoved(int position);
    void    on_obvsl_Smax_sliderMoved(int position);
    void    on_obvsl_Vmax_sliderMoved(int position);
    void    on_obhsl_sb_sliderMoved(int position);
    void    on_obhsl_vt_sliderMoved(int position);

    void    on_obhsl_sb_sliderPressed();
    void    on_obhsl_sb_sliderReleased();
    void    on_obhsl_vt_sliderPressed();
    void    on_obhsl_vt_sliderReleased();
    void    on_obvsl_Hmin_sliderPressed();
    void    on_obvsl_Hmin_sliderReleased();
    void    on_obvsl_Smin_sliderPressed();
    void    on_obvsl_Vmin_sliderReleased();
    void    on_obvsl_Smin_sliderReleased();
    void    on_obvsl_Vmin_sliderPressed();
    void    on_obvsl_Hmax_sliderPressed();
    void    on_obvsl_Hmax_sliderReleased();
    void    on_obvsl_Smax_sliderPressed();
    void    on_obvsl_Smax_sliderReleased();
    void    on_obvsl_Vmax_sliderPressed();
    void    on_obvsl_Vmax_sliderReleased();
    void    on_obhsl_frame_sliderPressed();
    void    on_obhsl_frame_sliderReleased();

    void    runSeries();
    void    on_runSeries();
    void    on_drawSeries();

    void    on_obckb_runseries_stateChanged(int arg1);

    double normalizeCoh(double dCoh, double dLbound);
    double coh4Adjacents(double dAngle, const cv::Point2f& pA, const cv::Point2f& pB, const cv::Point2f& pC, double dVt);

protected:
    QString qs_PathName;
    QTimer  qt_playVDO;
    cv::Mat m_frame, m_prevframe, m_framesequence, m_vector;
    int     i_frame, i_frameTotal;
    int     i_Hmin, i_Smin, i_Vmin;
    int     i_Hmax, i_Smax, i_Vmax;
    int     i_sb;
    double  d_vt;
    int     i_vDotR, i_vDotG, i_vDotB;
    int     i_vLineR, i_vLineG, i_vLineB;

    QTimer*     qt_runSeries;

    /// on_runSeries()
    int         i_begin;
    int         i_end;
    bool        b_firstReleased;

    double      d_StdRad;
    ///double      d_StdAng;
    double      d_StdCoh;
    ///double      d_StdMCoh;
    double      d_AvgRadius;
    ///double      d_AvgAngle;
    double      d_AvgCoh;
    ///double      d_AvgMCoh;

    double      md_Rad[10][30];
    double      md_Coh[10][30];
    double      md_covRad[10][30];
    double      md_covCoh[10][30];

    ///double      md_Ang[10][30];
    ///double      md_covAng[10][30];

    int         i_row;
    int         i_col;
    bool        b_out;

    int         i_prevVecCount;

    double      dAr_Rad[100];
    ///double      dAr_Ang[100];
    double      dAr_Coh[100];
    ///double      dAr_MCoh[100];

    ///int iAr_RadMag[21];
    ///int iAr_CohMag[101];

};
#endif // MAINWINDOW_H
