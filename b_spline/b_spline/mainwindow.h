#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "matrix.h"

#include <QKeyEvent>
#include <QMouseEvent>
#include <QIODevice>

#include <QDebug>
#include <QMessageBox>
#include <QPalette>
#include <qmath.h>
#include <QPainter>
#include <QImage>
#include <QPixmap>
#include <QBitmap>
#include <QPointF>
#include <QColor>
#include <QGraphicsScene>
#include <QGraphicsView>

#define FIT_POINT_NUM_LIMIT 300
#define CTRL_POINT_NUM_LIMIT (FIT_POINT_NUM_LIMIT + 2)
#define KNOT_NUM_LIMIT (CTRL_POINT_NUM_LIMIT + 4) // Knot_num = ctrl_num + level + 1
typedef double REAL;

typedef struct
{
    REAL X;
    REAL Y;
} COORDINATE_t;

typedef struct
{
    REAL X;
    REAL Y;
    REAL Width;
    REAL Height;
} RECT_MSG_t;

typedef struct
{
    RECT_MSG_t minRect;
    COORDINATE_t leftPos;
    COORDINATE_t rightPos;
    COORDINATE_t topPos;
    COORDINATE_t bottomPos;
    COORDINATE_t midPos;
} MIN_AREA_T;



typedef struct
{
    REAL value;
    REAL t;
} SPLINE_VALUE_t;

typedef struct
{
    MIN_AREA_T rect; //最小矩形
    REAL ratio;      //坐标分辨率
    REAL len;        //长度
    REAL maxSpeed;   //最大速度
    REAL bestTStep;  //最佳显示分辨率
    uint16_t fitNum; //拟合点数量
} SPLINE_MSG_T;

typedef struct
{
    REAL X;
    REAL Y;
    REAL crown;
} LWPOLYLINE_POINT_T;

typedef struct
{
    uint8_t isClose;
    uint8_t ctrlNum;
    LWPOLYLINE_POINT_T *ctrlVector;
} LWPOLYLINE_T;



#define BEZIER_SHOW_CTRL 1   //贝塞尔曲线控制点显示
#define B_SPLINE_SHOW_CTRL 0 // Bspline曲线控制点显示

#define _IN
#define _OUT

typedef struct MNode *PtrToMNode;
struct MNode
{
    int row;    //行
    int column; //列
    REAL **data;
};
typedef PtrToMNode Matrix;

typedef enum
{
    MATRIX_NO_ERROR = 0,                   //无错误
    MATRIX_INPUT_NO_SPECIFICATION,         //输入参数不规范
    MATRIX_FAILED_TO_ALLOCATE_HEAP_MEMORY, //分配堆内存失败
    MATRIX_ROWS_OR_COLUMNS_NOT_EQUAL,      //矩阵行数或列数不相等
    MATRIX_MULTIPLICATION,                 //矩阵乘法错误(第一个矩阵的列数不等于第二个矩阵行数)
    MATRIX_MUST_BE_SQUARE,                 //矩阵必须为方阵
    MATRIX_NOT_INVERTIBLE,                 //矩阵不可逆
    MATRIX_NOT_GAUSSORDINAL,               //矩阵不可高斯分解
} MATRIX_RES_T;

extern void Matrix_Zero(_IN Matrix mat);
extern Matrix Matrix_CreateZero(_IN int row, _IN int col);
extern void Matrix_Free(_IN Matrix mat);
extern Matrix Matrix_CreateEye(_IN int n);
extern void Matrix_Show(_IN Matrix mat);
extern void Matrix_SetData(_IN Matrix mat, _IN REAL data[]);
extern REAL Matrix_At(_IN Matrix mat, _IN int row, _IN int col);
extern void Matrix_Swap(_IN Matrix mat, _IN int n, _IN int m);
extern Matrix Matrix_Copy(_IN Matrix mat);
extern Matrix Matrix_Adjoint(_IN Matrix mat);
extern Matrix Matrix_Add(_IN Matrix mat_1, _IN Matrix mat_2);
extern Matrix Matrix_Sub(_IN Matrix mat_1, _IN Matrix mat_2);
extern Matrix Matrix_Transpose(_IN Matrix mat);
extern Matrix Matrix_Mult(_IN Matrix mat_1, _IN Matrix mat_2);
extern Matrix Matrix_MultConst(_IN Matrix mat, _IN REAL ratio);
extern REAL Matrix_Det(_IN Matrix mat);
extern Matrix Matrix_Inverse_LU(_IN Matrix mat);
extern Matrix Matrix_Inverse_EleTrans(_IN Matrix mat);

QT_BEGIN_NAMESPACE
namespace Ui
{
    class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    void Draw_Arc(int Xpos, int Ypos, int Radius, int startAngle, int stopAngle);
//    MIN_AREA_T GetMinArea_LwPolyLine(LWPOLYLINE_T *lwPolyLineMsg);
    SPLINE_VALUE_t Draw_BezierLine_GetMinAccel(uint8_t level, REAL tStep, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum);
    void Draw_BezierLine_FirstDerivative(uint8_t level, REAL tStep, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum);
    SPLINE_VALUE_t Draw_BezierLine_GetMaxSpeed(uint8_t level, REAL tStep, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum);
    MIN_AREA_T Draw_BezierLine(uint8_t level, REAL tStep, REAL ratio, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum);
    SPLINE_VALUE_t Draw_Bspline_GetMaxSpeedWithKnot(uint8_t level, REAL tStep, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum, REAL *KnotVector, uint16_t KnotNum);
    SPLINE_VALUE_t Draw_Bspline_GetMaxSpeedUneven(REAL tStep, uint8_t level, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum);
    SPLINE_VALUE_t Draw_Bspline_GetMaxSpeedUniform(REAL tStep, uint8_t level, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum);
    REAL Draw_Bspline_GetBidegN(uint16_t level, REAL *KnotVector, uint16_t KnotNum, uint16_t KnotIndex, REAL t);
    REAL Draw_Bspline_GetBidegP(uint16_t level, REAL *KnotVector, uint16_t KnotNum, uint16_t KnotIndex, REAL t);
    COORDINATE_t Draw_Bspline_DeBoor(uint16_t level,  REAL t, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum,REAL *KnotVector, uint16_t KnotNum);
    COORDINATE_t Draw_Bspline_DeBoor_FirstDerivative(uint16_t level,  REAL t,COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum,REAL *KnotVector, uint16_t KnotNum);
    uint16_t Draw_Bspline_CreateKnotVectorClamped(uint16_t level, uint16_t ctrlPointNum, REAL *KnotVector);
    uint16_t Draw_Bspline_CreateKnotVectorSumChord(uint16_t level, COORDINATE_t *fitPointVector, uint16_t fitPointNum, REAL *KnotVector);
    MIN_AREA_T Draw_BsplineWithKnot(uint16_t level, REAL tStep, REAL ratio, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum, REAL *KnotVector, uint16_t KnotNum);
    MIN_AREA_T Draw_BsplineClamped(uint16_t level, REAL tStep, REAL ratio, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum);
    MIN_AREA_T Draw_BsplineWithFit(uint16_t level, REAL ratio, COORDINATE_t *FitVector, uint16_t FitNum);
    uint16_t Draw_Bspline_GetFitVector(uint16_t level,REAL ratio, REAL speed, COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum, REAL *KnotVector, uint16_t KnotNum, COORDINATE_t *fitPointVector);
    ~MainWindow();

private slots:
    void on_pushButton_clear_clicked();
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveOnCurve(QMouseEvent *event);

private:
    Ui::MainWindow *ui;
    QList<QString> mPatterList;
    QGraphicsScene *my_scene;
    QPixmap image;
    int ctrlPointCnt = 0;
    QList<QPoint> ctrlPonitList;
    COORDINATE_t ctrlPointForBspline[300];
};
#endif // MAINWINDOW_H
