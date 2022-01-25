#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->canvas, SIGNAL(mouseMove(QMouseEvent *)), this, SLOT(mouseMoveOnCurve(QMouseEvent *)));     //! [鼠标移动事件]
    connect(ui->canvas, SIGNAL(mousePress(QMouseEvent *)), this, SLOT(mousePressEvent(QMouseEvent *)));     //! [鼠标按压事件]
    connect(ui->canvas, SIGNAL(mouseRelease(QMouseEvent *)), this, SLOT(mouseReleaseEvent(QMouseEvent *))); //! [鼠标释放事件]

    my_scene = new QGraphicsScene; //新建绘图选项
    image = QPixmap(800, 480);
    image.fill(qRgb(255, 255, 255)); //将位图背景设置为白色
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
    mPatterList << "多段三阶贝塞尔" << QString::fromLocal8Bit("B-spline-ctrl") << QString::fromLocal8Bit("B-spline-fit") << QString::fromLocal8Bit("B-spline-get-fit");
    ui->comboBox_lineType->addItems(mPatterList);
    ui->comboBox_lineType->setCurrentIndex(0);
    memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));

//    LWPOLYLINE_POINT_T ctrlPoint[30] = {
//        {199, 50, 0},
//        {102.45, 50, 1},
//        {100, 53.45, 0},
//        {100, 97.56, 0},
//        {102.55, 100, 0},
//        {197.55, 100, -0.41421},
//        {200, 97.56, 0},
//        {200, 51, -0.41421},
//    };
//    LWPOLYLINE_T lwMsg = {1, 8, ctrlPoint};
//    GetMinArea_LwPolyLine(&lwMsg);
}

MainWindow::~MainWindow()
{
    delete ui;
}

/**
 * @brief 初始化矩阵(将所有元素初始化为０)
 * @param {_IN Matrix} mat
 * @return {*}
 * @note
 */
void Matrix_Zero(_IN Matrix mat)
{
    int i, j;
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
        {
            mat->data[i][j] = 0;
        }
    }
}

/**
 * @brief 　创建0矩阵
 * @param {_IN int} row
 * @param {_IN int} col
 * @return {*}
 * @note
 */
Matrix Matrix_CreateZero(_IN int row, _IN int col)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    if (row <= 0 || col <= 0)
    {
        res = MATRIX_INPUT_NO_SPECIFICATION;
        return NULL;
    }
    Matrix mat = (Matrix)malloc(sizeof(struct MNode)); //　分配结构体指针
    if (row > 0 && col > 0)
    {
        mat->row = row;
        mat->column = col;
        mat->data = (REAL **)malloc(row * sizeof(REAL *)); // 分配头指针
        if (mat->data == NULL)
        {
            res = MATRIX_FAILED_TO_ALLOCATE_HEAP_MEMORY;
            return NULL;
        }
        int i;
        for (i = 0; i < row; i++)
        {
            *(mat->data + i) = (REAL *)malloc(col * sizeof(REAL)); //　分配每行的指针
            if (mat->data[i] == NULL)
            {
                res = MATRIX_FAILED_TO_ALLOCATE_HEAP_MEMORY;
                return NULL;
            }
        }
        Matrix_Zero(mat);
    }
    return mat;
}

/**
 * @brief 释放申请的矩阵空间
 * @param {_IN Matrix} mat
 * @return {*}
 * @note
 */
void Matrix_Free(_IN Matrix mat)
{
    for (int i = 0; i < mat->row; i++) // 释放行指针
    {
        free(mat->data[i]);
    }
    free(mat->data); // 释放头指针
    free(mat);       // 释放结构体指针
}

/**
 * @brief 　创建单位矩阵
 * @param {int} n
 * @return {*}
 * @note
 */
Matrix Matrix_CreateEye(_IN int n)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    int i, j;
    Matrix mat;
    if (n <= 0)
    {
        res = MATRIX_INPUT_NO_SPECIFICATION;
        return NULL;
    }
    mat = Matrix_CreateZero(n, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
                mat->data[i][j] = 1;
            else
                mat->data[i][j] = 0;
        }
    }
    return mat;
}

/**
 * @brief 　打印矩阵
 * @param {_IN Matrix} mat
 * @return {*}
 * @note
 */
void Matrix_Show(_IN Matrix mat)
{
    int i, j;
    //    printf("%s\n", s);
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
        {
            printf("%.6f\t", mat->data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * @brief 给矩阵每个元素赋值
 * @param {_IN Matrix} mat
 * @param {_IN REAL} data
 * @return {*}
 * @note
 */
void Matrix_SetData(_IN Matrix mat, _IN REAL data[])
{
    int i, j;
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
        {
            mat->data[i][j] = data[i * mat->column + j];
        }
    }
}

/**
 * @brief 　取出矩阵某行某列的元素
 * @param {_IN Matrix} mat
 * @param {_IN int} row
 * @param {_IN int} col
 * @return {*}
 * @note
 */
REAL Matrix_At(_IN Matrix mat, _IN int row, _IN int col)
{
    REAL rst;
    if (row <= 0 || col <= 0)
    {
        return 0;
    }
    rst = mat->data[row - 1][col - 1];
    return rst;
}

/**
 * @brief 矩阵第n行与第m行互换
 * @param {_IN Matrix} mat
 * @param {_IN int} n
 * @param {_IN int} m
 * @return {*}
 * @note
 */
void Matrix_Swap(_IN Matrix mat, _IN int n, _IN int m)
{
    REAL temp;
    for (int i = 0; i < mat->column; i++)
    {
        temp = mat->data[n - 1][i];
        mat->data[n - 1][i] = mat->data[m - 1][i];
        mat->data[m - 1][i] = temp;
    }
}

/**
 * @brief 　对一个矩阵进行复制
 * @param {_IN Matrix} mat
 * @return {*}
 * @note
 */
Matrix Matrix_Copy(_IN Matrix mat)
{
    Matrix matCpy;
    matCpy = Matrix_CreateZero(mat->row, mat->column);
    for (int i = 0; i < mat->row; i++)
    {
        for (int j = 0; j < mat->column; j++)
            matCpy->data[i][j] = mat->data[i][j];
    }
    return matCpy;
}

/**
 * @brief 伴随矩阵
 * @param {_IN Matrix} mat
 * @return {*}
 * @note
 */
Matrix Matrix_Adjoint(_IN Matrix mat)
{
    //    Matrix adj;
    return mat;
}

/**
 * @brief 矩阵加法
 * @param {_IN Matrix} mat_1
 * @param {_IN Matrix} mat_2
 * @param {_OUT Matrix} matRes
 * @return {*}
 * @note
 */
Matrix Matrix_Add(_IN Matrix mat_1, _IN Matrix mat_2)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    Matrix matRes;
    if ((mat_1->column != mat_2->column) || (mat_1->row != mat_2->row))
    {
        res = MATRIX_ROWS_OR_COLUMNS_NOT_EQUAL;
        return NULL;
    }
    matRes = Matrix_CreateZero(mat_1->row, mat_1->column);
    for (int i = 0; i < mat_1->row; i++)
    {
        for (int j = 0; j < mat_1->column; j++)
            matRes->data[i][j] = mat_1->data[i][j] + mat_2->data[i][j];
    }
    return matRes;
}

/**
 * @brief 矩阵减法
 * @param {_IN Matrix} mat_1
 * @param {_IN Matrix} mat_2
 * @return {*}
 * @note
 */
Matrix Matrix_Sub(_IN Matrix mat_1, _IN Matrix mat_2)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    Matrix matRes;
    if ((mat_1->column != mat_2->column) || (mat_1->row != mat_2->row))
    {
        res = MATRIX_ROWS_OR_COLUMNS_NOT_EQUAL;
        return NULL;
    }
    matRes = Matrix_CreateZero(mat_1->row, mat_1->column);
    for (int i = 0; i < mat_1->row; i++)
    {
        for (int j = 0; j < mat_1->column; j++)
            matRes->data[i][j] = mat_1->data[i][j] - mat_2->data[i][j];
    }
    return matRes;
}

/**
 * @brief 矩阵转置
 * @param {_IN Matrix} mat
 * @param {_OUT Matrix} matRes
 * @return {*}
 * @note
 */
Matrix Matrix_Transpose(_IN Matrix mat)
{
    Matrix matRes;
    matRes = Matrix_CreateZero(mat->row, mat->column);
    for (int i = 0; i < mat->row; i++)
    {
        for (int j = 0; j < mat->column; j++)
            matRes->data[i][j] = mat->data[i][j];
    }
    return matRes;
}

/**
 * @brief 矩阵乘法
 * @param {_IN Matrix} mat_1
 * @param {_IN Matrix} mat_2
 * @return {*}
 * @note
 */
Matrix Matrix_Mult(_IN Matrix mat_1, _IN Matrix mat_2)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    Matrix matRes;
    if (mat_1->column != mat_2->row)
    {
        res = MATRIX_MULTIPLICATION;
        return NULL;
    }
    else
    {
        matRes = Matrix_CreateZero(mat_1->row, mat_2->column);
        for (int i = 0; i < mat_1->row; i++)
        {
            for (int j = 0; j < mat_2->column; j++)
            {
                for (int m = 0; m < mat_1->column; m++)
                    matRes->data[i][j] += mat_1->data[i][m] * mat_2->data[m][j];
            }
        }
    }
    return matRes;
}

/**
 * @brief 矩阵乘法常数
 * @param {_IN Matrix} mat
 * @param {_IN REAL} ratio
 * @return {*}
 * @note
 */
Matrix Matrix_MultConst(_IN Matrix mat, _IN REAL ratio)
{
    for (int i = 0; i < mat->row; i++)
    {
        for (int j = 0; j < mat->column; j++)
        {
            mat->data[i][j] *= ratio;
        }
    }
    return mat;
}

/**
 * @brief 通过行变换（列主元消去法）将矩阵变换成上三角矩阵
 * @param {_IN Matrix} mat
 * @return {*}
 * @note 注意行变换会改变行列式的符号
 */
REAL Matrix_Det(_IN Matrix mat)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    REAL det = 1.0;
    Matrix mat_;
    if (mat->row != mat->column)
    {
        res = MATRIX_MUST_BE_SQUARE;
        return NULL;
    }
    if (mat->row == 1 && mat->column == 1)
    {
        return mat->data[0][0];
    }
    // 为防止改变原矩阵，此处进行复制处理
    mat_ = Matrix_Copy(mat);
    int n = mat_->row, s = 0, cnt = 0;
    REAL L, a_max;
    for (int k = 0; k < n - 1; k++)
    {
        cnt = k;
        a_max = fabs(mat_->data[k][k]);
        for (int i = k; i < n; i++)
        {
            if (fabs(mat_->data[i][k]) > a_max)
            {
                a_max = fabs(mat_->data[i][k]);
                cnt = i; // 确定下标i
            }
        }
        //　换行
        REAL temp;
        if (cnt != k)
        {
            s++; // 换行次数
            for (int j = 0; j < n; j++)
            {
                temp = mat_->data[cnt][j];
                mat_->data[cnt][j] = mat_->data[k][j];
                mat_->data[k][j] = temp;
            }
        }
        // 消元计算
        for (int i = k + 1; i < n; i++)
        {
            L = mat_->data[i][k] / mat_->data[k][k];
            for (int j = k + 1; j < n; j++)
                mat_->data[i][j] = mat_->data[i][j] - L * mat_->data[k][j];
        }
    }
    if (s % 2 == 0)
        s = 1;
    else
        s = -1;
    for (int i = 0; i < n; i++)
        det *= mat_->data[i][i];
    det *= s;
    Matrix_Free(mat_); //释放掉复制矩阵的内存
    return det;
}

/**
 * @brief 矩阵求逆，利用矩阵LU分解求逆（直接三角分解）
 * @param {_IN Matrix} mat
 * @return {*}
 * @note
 */
Matrix Matrix_Inverse_LU(_IN Matrix mat)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    Matrix L, U;
    Matrix L_, U_;
    Matrix matRes;
    if (mat->column != mat->row)
    {
        res = MATRIX_MUST_BE_SQUARE;
        return NULL;
    }
    if (Matrix_Det(mat) == 0)
    {
        res = MATRIX_NOT_INVERTIBLE;
        return NULL;
    }
    // Ｌ矩阵和Ｕ矩阵

    int n = mat->column;
    matRes = Matrix_CreateZero(n, n);
    L = Matrix_CreateZero(n, n);
    U = Matrix_CreateZero(n, n);
    // 计算Ｕ的第一行元素
    for (int j = 0; j < n; j++)
        U->data[0][j] = mat->data[0][j];
    //　计算Ｌ的第一列元素
    for (int i = 0; i < n; i++)
    {
        L->data[i][0] = mat->data[i][0] / U->data[0][0];
        L->data[i][i] = 1.0;
    }
    REAL sum_u = 0, sum_l = 0;
    //　分别计算Ｕ和Ｌ的第２到ｎ行、列元素
    for (int k = 1; k < n; k++)
    {
        // 求U，U矩阵从第２行迭代到第n行，且U矩阵先于L矩阵一个节拍
        for (int j = k; j < n; j++)
        {
            sum_u = 0;
            for (int m = 0; m <= k - 1; m++)
                sum_u += L->data[k][m] * U->data[m][j];
            U->data[k][j] = mat->data[k][j] - sum_u;
        }
        //　求Ｌ，L矩阵从第２列迭代到第n列
        for (int i = k + 1; i < n; i++)
        {
            sum_l = 0;
            for (int m = 0; m <= k - 1; m++)
                sum_l += L->data[i][m] * U->data[m][k];
            L->data[i][k] = (mat->data[i][k] - sum_l) / U->data[k][k];
        }
    }
    // 分别求下三角矩阵Ｌ和上三角矩阵Ｕ的逆矩阵

    L_ = Matrix_CreateZero(n, n);
    U_ = Matrix_CreateZero(n, n);
    // 求矩阵Ｕ的逆
    REAL sum_u_;
    for (int i = 0; i < n; i++)
    {
        U_->data[i][i] = 1.0 / U->data[i][i]; // 对角线元素的值直接取倒数
        for (int j = i - 1; j >= 0; j--)
        {
            sum_u_ = 0;
            for (int k = j + 1; k <= i; k++)
                sum_u_ += U->data[j][k] * U_->data[k][i];
            U_->data[j][i] = -sum_u_ / U->data[j][j]; // 迭代计算，按列倒序依次得到每一个值
        }
    }
    // 求Ｌ的逆
    for (int i = 0; i < n; i++)
    {
        L_->data[i][i] = 1; // 对角线元素的值直接取倒数，这里为１
        for (int k = i + 1; k < n; k++)
        {
            for (int j = i; j <= k - 1; j++)
                L_->data[k][i] -= L->data[k][j] * L_->data[j][i]; // 迭代计算，按列顺序依次得到每一个值
        }
    }
    // 已经得到Ｌ和Ｕ的逆矩阵
    matRes = Matrix_Mult(U_, L_);

    Matrix_Free(L_);
    Matrix_Free(U_);
    Matrix_Free(L);
    Matrix_Free(U);
    return matRes;
}

/**
 * @brief 矩阵求逆，利用初等行变换求逆
 * @param {_IN Matrix} mat
 * @return {*}
 * @note
 */
Matrix Matrix_Inverse_EleTrans(_IN Matrix mat)
{
    MATRIX_RES_T res = MATRIX_NO_ERROR;
    Matrix matRes;
    Matrix mat_;
    if (mat->column != mat->row)
    {
        res = MATRIX_MUST_BE_SQUARE;
        return NULL;
    }
    if (Matrix_Det(mat) == 0)
    {
        res = MATRIX_NOT_INVERTIBLE;
        return NULL;
    }
    int n = mat->row;
    // 为防止改变原矩阵，此处进行复制处理
    mat_ = Matrix_Copy(mat);
    //　创建单位矩阵
    matRes = Matrix_CreateEye(n);
    REAL e, a_max;
    int i, j, k, t, cnt;
    for (k = 0; k < n - 1; k++)
    {
        a_max = fabs(mat_->data[k][k]);
        cnt = k;
        // 选主元
        for (i = k; i < n; i++)
        {
            if (fabs(mat_->data[i][k]) > a_max)
            {
                a_max = fabs(mat_->data[i][k]);
                cnt = i;
            }
        }
        //　换行，原矩阵换行的同时，单位矩阵也换行
        REAL temp, temp_e;
        if (cnt != k)
        {
            for (j = 0; j < n; j++)
            {
                temp = mat_->data[k][j];
                mat_->data[k][j] = mat_->data[cnt][j];
                mat_->data[cnt][j] = temp;
                temp_e = matRes->data[k][j];
                matRes->data[k][j] = matRes->data[cnt][j];
                matRes->data[cnt][j] = temp_e;
            }
        }
        // 消元
        for (i = k + 1; i < n; i++)
        {
            e = mat_->data[i][k] / mat_->data[k][k];
            for (j = 0; j < n; j++)
            {
                mat_->data[i][j] = mat_->data[i][j] - e * mat_->data[k][j];
                matRes->data[i][j] = matRes->data[i][j] - e * matRes->data[k][j];
            }
        }
    }
    // 将矩阵每行的行首元素化为１
    for (i = 0; i < n; i++)
    {
        e = 1.0 / mat_->data[i][i];
        for (j = 0; j < n; j++)
        {
            mat_->data[i][j] = mat_->data[i][j] * e;
            matRes->data[i][j] = matRes->data[i][j] * e;
        }
    }
    // 从最后一排开始消元，把增广矩阵左边的矩阵化为单位矩阵
    for (i = n - 1; i > 0; i--)
    {
        for (j = i - 1; j >= 0; j--)
        {
            e = mat_->data[j][i] / mat_->data[i][i];
            for (t = 0; t < n; t++)
            {
                mat_->data[j][t] = mat_->data[j][t] - e * mat_->data[i][t];
                matRes->data[j][t] = matRes->data[j][t] - e * matRes->data[i][t];
            }
        }
    }
    Matrix_Free(mat_);
    return matRes;
}

//! [鼠标按压事件]
void MainWindow::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton)
    {
        if (ui->comboBox_lineType->currentIndex() == 0)
        {
            if (ctrlPointCnt < 4)
            {
                QPoint p_ab = event->globalPos();
                QPoint p_re = event->pos();
                if (p_re.x() < ui->canvas->pos().x() ||
                        p_re.x() > ui->canvas->pos().x() + ui->canvas->width() ||
                        p_re.y() < ui->canvas->pos().y() ||
                        p_re.y() > ui->canvas->pos().y() + ui->canvas->height())
                {
                    return;
                }
                p_re.setX(p_re.x() - ui->canvas->pos().x() - 15);
                p_re.setY(p_re.y() - ui->canvas->pos().y() - 15);
                ctrlPointCnt++;

                QPainter painter(&image); //选入绘图设备中。
                painter.setPen(Qt::red);

                QRectF rectangle(p_re.x(), p_re.y(), 3, 3);
                painter.drawEllipse(rectangle);
                painter.end();
                my_scene->addPixmap(image);
                ui->canvas->setScene(my_scene);
                ui->canvas->show();
                QString msg;
                msg += "点" + QString::number(ctrlPointCnt) + ":(" + QString::number(p_re.x(), 'g', 6) +
                        "," + QString::number(p_re.y(), 'g', 6) + ")\n";
                ui->textBrowser_msg->append(msg);
                ui->textBrowser_msg->show();
                ui->comboBox_lineType->setEnabled(false);
                ctrlPonitList.append(p_re);
                if (ctrlPointCnt == 4)
                {
                    COORDINATE_t pointList[4];
                    pointList[0].X = ctrlPonitList.at(0).x();
                    pointList[0].Y = ctrlPonitList.at(0).y();
                    pointList[1].X = ctrlPonitList.at(1).x();
                    pointList[1].Y = ctrlPonitList.at(1).y();
                    pointList[2].X = ctrlPonitList.at(2).x();
                    pointList[2].Y = ctrlPonitList.at(2).y();
                    pointList[3].X = ctrlPonitList.at(3).x();
                    pointList[3].Y = ctrlPonitList.at(3).y();
                    Draw_BezierLine(3, 0, 1, pointList, 4);
                    ctrlPonitList.clear();
                    ui->textBrowser_msg->clear();
                    ui->textBrowser_msg->show();
                    ctrlPointCnt = 0;
                    ui->comboBox_lineType->setEnabled(true);
                    memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));
                }
            }
            else
            {

                ui->textBrowser_msg->clear();
                ui->textBrowser_msg->show();
                ctrlPointCnt = 0;
                ui->comboBox_lineType->setEnabled(true);
            }
        }
        else
        {
            if (ctrlPointCnt < CTRL_POINT_NUM_LIMIT)
            {
                QPoint p_ab = event->globalPos();
                QPoint p_re = event->pos();
                if (p_re.x() < ui->canvas->pos().x() ||
                        p_re.x() > ui->canvas->pos().x() + ui->canvas->width() ||
                        p_re.y() < ui->canvas->pos().y() ||
                        p_re.y() > ui->canvas->pos().y() + ui->canvas->height())
                {
                    return;
                }
                p_re.setX(p_re.x() - ui->canvas->pos().x() - 15);
                p_re.setY(p_re.y() - ui->canvas->pos().y() - 15);

                QPainter painter(&image); //选入绘图设备中。
                painter.setPen(Qt::red);

                QRectF rectangle(p_re.x(), p_re.y(), 3, 3);
                painter.drawEllipse(rectangle);
                painter.end();
                my_scene->addPixmap(image);
                ui->canvas->setScene(my_scene);
                ui->canvas->show();
                QString msg;
                msg += "点" + QString::number(ctrlPointCnt) + ":(" + QString::number(p_re.x(), 'g', 6) +
                        "," + QString::number(p_re.y(), 'g', 6) + ")\n";
                ui->textBrowser_msg->append(msg);
                ui->textBrowser_msg->show();

                ctrlPointForBspline[ctrlPointCnt].X = p_re.x();
                ctrlPointForBspline[ctrlPointCnt].Y = p_re.y();
                ctrlPointCnt++;
                ui->comboBox_lineType->setEnabled(false);
            }
        }
    }
    else if (event->button() == Qt::RightButton)
    {
        if (ui->comboBox_lineType->currentIndex() == 1)
        {
            if (ctrlPointCnt >= 4)
            {
                Draw_BsplineClamped(3, 0, 1, ctrlPointForBspline, ctrlPointCnt);
                ctrlPointCnt = 0;
                ui->comboBox_lineType->setEnabled(true);
                memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));
            }
            else if (ctrlPointCnt == 3)
            {
                Draw_BsplineClamped(2, 0, 1, ctrlPointForBspline, ctrlPointCnt);
                ctrlPointCnt = 0;
                ui->comboBox_lineType->setEnabled(true);
                memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));
            }
        }
        else if (ui->comboBox_lineType->currentIndex() == 2)
        {
            if (ctrlPointCnt >= 2)
            {
                Draw_BsplineWithFit(3, 1, ctrlPointForBspline, ctrlPointCnt);
                ctrlPointCnt = 0;
                ui->comboBox_lineType->setEnabled(true);
                memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));
            }
        }
        else if (ui->comboBox_lineType->currentIndex() == 3)
        {
            if (ctrlPointCnt >= 4)
            {
                REAL KnotVector[CTRL_POINT_NUM_LIMIT];
                COORDINATE_t fitVector[CTRL_POINT_NUM_LIMIT];
                uint16_t KnotNum;
                uint16_t msg;
                KnotNum = Draw_Bspline_CreateKnotVectorClamped(3, ctrlPointCnt, KnotVector);
                msg = Draw_Bspline_GetFitVector(3, 1, 50, ctrlPointForBspline, ctrlPointCnt, KnotVector, KnotNum, fitVector);

                // Draw_BsplineWithFit(3,1,fitVector,msg);
                ctrlPointCnt = 0;
                ui->comboBox_lineType->setEnabled(true);
                memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));
            }
            else if (ctrlPointCnt == 3)
            {
                REAL KnotVector[CTRL_POINT_NUM_LIMIT];
                COORDINATE_t fitVector[CTRL_POINT_NUM_LIMIT];
                uint16_t KnotNum;
                uint16_t msg;
                KnotNum = Draw_Bspline_CreateKnotVectorClamped(2, ctrlPointCnt, KnotVector);
                msg = Draw_Bspline_GetFitVector(2, 1, 50, ctrlPointForBspline, ctrlPointCnt, KnotVector, KnotNum, fitVector);

                // Draw_BsplineWithFit(3,1,fitVector,msg);
                ctrlPointCnt = 0;
                ui->comboBox_lineType->setEnabled(true);
                memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));
            }
        }
    }
}

//! [鼠标释放事件]
void MainWindow::mouseReleaseEvent(QMouseEvent *event)
{
}

//! [鼠标移动事件]
void MainWindow::mouseMoveOnCurve(QMouseEvent *event)
{
    QPoint p_ab = event->globalPos();
    QPoint p_re = event->pos();
    ui->textBrowser_X->clear();
    ui->textBrowser_Y->clear();
    ui->textBrowser_X->append(QString::number(p_re.x(), 'g', 6));
    ui->textBrowser_Y->append(QString::number(p_re.y(), 'g', 6));
    ui->textBrowser_X->show();
    ui->textBrowser_Y->show();
}

int ftoint(REAL a)
{
    return (int)(a + 0.5);
}

/**
 * @brief 绘制圆弧
 * @param {float} Xpos
 * @param {float} Ypos
 * @param {float} Radius
 * @param {float} startAngle
 * @param {float} stopAngle
 * @return {*}
 * @note
 */
void MainWindow::Draw_Arc(int Xpos, int Ypos, int Radius, int startAngle, int stopAngle)
{
    int x = -Radius, y = 0, err = 2 - 2 * Radius, e2;
    uint8_t dirArc = 0; //左方
    QPainter painter(&image);        //选入绘图设备中。
    painter.setPen(Qt::blue);
    if (startAngle == stopAngle)
    {
        //        LCD_DrawCircle(Xpos, Ypos, Radius);
        return;
    }
    while (startAngle > 360)
    {
        startAngle -= 360;
    }
    while (startAngle < 0)
    {
        startAngle += 360;
    }
    while (stopAngle > 360)
    {
        stopAngle -= 360;
    }
    while (stopAngle < 0)
    {
        stopAngle += 360;
    }
    if (startAngle > stopAngle)
    {
        dirArc = 1; //右方
    }

    do
    {

        int tmpX = x;
        int tmpY = y;
        int tmpArccos = 0;
        tmpArccos = (int)(acos((float)tmpX / (float)Radius) * 180 / M_PI);
        if (tmpY < 0)
        {
            tmpArccos = 360 - tmpArccos;
        }
        if (dirArc == 1)
        {
            if (tmpArccos <= stopAngle || tmpArccos >= startAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }
        else
        {
            if (tmpArccos >= startAngle && tmpArccos <= stopAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }

        tmpX = -x;
        tmpY = y;
        tmpArccos = 0;
        tmpArccos = (int)(acos((float)tmpX / (float)Radius) * 180 / M_PI);
        if (tmpY < 0)
        {
            tmpArccos = 360 - tmpArccos;
        }
        if (dirArc == 1)
        {
            if (tmpArccos <= stopAngle || tmpArccos >= startAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }
        else
        {
            if (tmpArccos >= startAngle && tmpArccos <= stopAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }
        tmpX = x;
        tmpY = -y;
        tmpArccos = 0;
        tmpArccos = (int)(acos((float)tmpX / (float)Radius) * 180 / M_PI);
        if (tmpY < 0)
        {
            tmpArccos = 360 - tmpArccos;
        }
        if (dirArc == 1)
        {
            if (tmpArccos <= stopAngle || tmpArccos >= startAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }
        else
        {
            if (tmpArccos >= startAngle && tmpArccos <= stopAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }
        tmpX = -x;
        tmpY = -y;
        tmpArccos = 0;
        tmpArccos = (int)(acos((float)tmpX / (float)Radius) * 180 / M_PI);
        if (tmpY < 0)
        {
            tmpArccos = 360 - tmpArccos;
        }
        if (dirArc == 1)
        {
            if (tmpArccos <= stopAngle || tmpArccos >= startAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }
        else
        {
            if (tmpArccos >= startAngle && tmpArccos <= stopAngle)
            {
                //                LCD_PutPixel((uint16_t)(Xpos + tmpX), (uint16_t)(Ypos + tmpY));
                painter.drawPoint(QPoint((Xpos + tmpX), (Ypos + tmpY)));
            }
            else
            {
            }
        }

        e2 = err;
        if (e2 <= y)
        {
            err += ++y * 2 + 1;
            if (-x == y && e2 <= x)
                e2 = 0;
        }
        if (e2 > x)
            err += ++x * 2 + 1;
    } while (x <= 0);

    painter.end();
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
}

//MIN_AREA_T MainWindow::GetMinArea_LwPolyLine(LWPOLYLINE_T *lwPolyLineMsg)
//{
//    MIN_AREA_T msg = {
//        .minRect = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y, 0, 0},
//        .leftPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .rightPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .topPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .midPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .bottomPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y}};
//    MIN_AREA_T msgTmp = {
//        .minRect = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y, 0, 0},
//        .leftPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .rightPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .topPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .midPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y},
//        .bottomPos = {lwPolyLineMsg->ctrlVector[0].X, lwPolyLineMsg->ctrlVector[0].Y}};

//    REAL Radius = 0;
//    COORDINATE_t org = {0, 0};
//    REAL angleA = 0, angleB = 0;
//    uint8_t direction = 0; // 0：正向 1：负向

//    REAL k = 0, _k = 0;
//    REAL _b = 0;
//    REAL h = 0;

//    REAL a = 0, b = 0, c = 0;
//    REAL delta = 0;
//    REAL x1 = 0, x2 = 0;

//    LWPOLYLINE_POINT_T PointNext = {0};
//    for (uint16_t i = 0; i < lwPolyLineMsg->ctrlNum; i++)
//    {

//        if (i == lwPolyLineMsg->ctrlNum - 1)
//        {
//            if (lwPolyLineMsg->isClose)
//            {
//                PointNext = lwPolyLineMsg->ctrlVector[0];
//            }
//            else
//            {
//                break;
//            }
//        }
//        else
//        {
//            PointNext = lwPolyLineMsg->ctrlVector[i + 1];
//        }
//        if (lwPolyLineMsg->ctrlVector[i].crown == 0)
//        {
////            msgTmp = GetMinArea_Line(lwPolyLineMsg->ctrlVector[i].X, lwPolyLineMsg->ctrlVector[i].Y,
////                                     PointNext.X, PointNext.Y);
//        }
//        else
//        {

//            /* 求解方向 */
//            if (lwPolyLineMsg->ctrlVector[i].X <= lwPolyLineMsg->ctrlVector[i + 1].X && lwPolyLineMsg->ctrlVector[i].Y < lwPolyLineMsg->ctrlVector[i + 1].Y)
//            {
//                direction = 0;
//            }
//            else if (lwPolyLineMsg->ctrlVector[i].X < lwPolyLineMsg->ctrlVector[i + 1].X && lwPolyLineMsg->ctrlVector[i].Y <= lwPolyLineMsg->ctrlVector[i + 1].Y)
//            {
//                direction = 0;
//            }
//            else
//            {
//                direction = 1;
//            }

//            /* 求解圆心坐标 */
//            /* 平行X轴 */
//            if (lwPolyLineMsg->ctrlVector[i].Y == PointNext.Y)
//            {
//                org.X = (lwPolyLineMsg->ctrlVector[i].X + PointNext.X) / 2;
//                org.Y = lwPolyLineMsg->ctrlVector[i].Y;
//                Radius = fabs(lwPolyLineMsg->ctrlVector[i].X - PointNext.X) / 2;
//                h = Radius * fabs(lwPolyLineMsg->ctrlVector[i].crown);
//                Radius = sqrt(h * h + Radius * Radius);
//                if (direction == 1)
//                {
//                    if (lwPolyLineMsg->ctrlVector[i].crown > 0)
//                    {
//                        org.X = org.X;
//                        org.Y += Radius;
//                    }
//                    else
//                    {
//                        org.X = org.X;
//                        org.Y -= Radius;
//                    }
//                }
//                else
//                {
//                    if (lwPolyLineMsg->ctrlVector[i].crown < 0)
//                    {
//                        org.X = org.X;
//                        org.Y += Radius;
//                    }
//                    else
//                    {
//                        org.X = org.X;
//                        org.Y -= Radius;
//                    }
//                }
//            }
//            /* 平行y轴 */
//            else if (lwPolyLineMsg->ctrlVector[i].X == PointNext.X)
//            {
//                org.X = lwPolyLineMsg->ctrlVector[i].X;
//                org.Y = (lwPolyLineMsg->ctrlVector[i].Y + PointNext.Y) / 2;
//                Radius = fabs(lwPolyLineMsg->ctrlVector[i].Y - PointNext.Y) / 2;
//                h = Radius * fabs(lwPolyLineMsg->ctrlVector[i].crown);
//                Radius = sqrt(h * h + Radius * Radius);
//                if (lwPolyLineMsg->ctrlVector[i].Y > PointNext.Y)
//                {
//                    if (lwPolyLineMsg->ctrlVector[i].crown > 0)
//                    {
//                        org.Y = org.Y;
//                        org.X += Radius;
//                    }
//                    else
//                    {
//                        org.Y = org.Y;
//                        org.X -= Radius;
//                    }
//                }
//                else
//                {
//                    if (lwPolyLineMsg->ctrlVector[i].crown < 0)
//                    {
//                        org.Y = org.Y;
//                        org.X += Radius;
//                    }
//                    else
//                    {
//                        org.Y = org.Y;
//                        org.X -= Radius;
//                    }
//                }
//            }
//            //任意圆弧
//            else
//            {
//                if (lwPolyLineMsg->ctrlVector[i].crown == 1 || lwPolyLineMsg->ctrlVector[i].crown == -1)
//                {

//                    Radius = sqrt((lwPolyLineMsg->ctrlVector[i].X - PointNext.X) *
//                                      (lwPolyLineMsg->ctrlVector[i].X - PointNext.X) +
//                                  (lwPolyLineMsg->ctrlVector[i].Y - PointNext.Y) *
//                                      (lwPolyLineMsg->ctrlVector[i].Y - PointNext.Y)) /
//                             2;
//                    org.X = (lwPolyLineMsg->ctrlVector[i].X + PointNext.X) / 2;
//                    org.Y = (lwPolyLineMsg->ctrlVector[i].Y + PointNext.Y) / 2;
//                }
//                else
//                {

//                    k = (lwPolyLineMsg->ctrlVector[i].Y - PointNext.Y) /
//                        (lwPolyLineMsg->ctrlVector[i].X - PointNext.X);
//                    _k = -1 / k;
//                    org.X = (lwPolyLineMsg->ctrlVector[i].X + PointNext.X) / 2;
//                    org.Y = (lwPolyLineMsg->ctrlVector[i].Y + PointNext.Y) / 2;
//                    _b = org.Y - _k * org.X;
//                    Radius = sqrt((lwPolyLineMsg->ctrlVector[i].X - PointNext.X) *
//                                      (lwPolyLineMsg->ctrlVector[i].X - PointNext.X) +
//                                  (lwPolyLineMsg->ctrlVector[i].Y - PointNext.Y) *
//                                      (lwPolyLineMsg->ctrlVector[i].Y - PointNext.Y));
//                    h = Radius / 2 * fabs(lwPolyLineMsg->ctrlVector[i].crown);
//                    Radius = (h * h + (Radius / 2) * (Radius / 2)) / (2 * h);
//                    if (fabs(lwPolyLineMsg->ctrlVector[i].crown) > 1)
//                    {
//                        h = h - Radius;
//                    }
//                    else
//                    {
//                        h = Radius - h;
//                    }

//                    a = (1 + _k * _k);
//                    b = (-2 * org.X + 2 * _k * _b - 2 * _k * org.Y);
//                    c = org.X * org.X + _b * _b - 2 * _b * org.Y + org.Y * org.Y - h * h;
//                    delta = b * b - 4 * a * c;
//                    if (delta < 0)
//                    {
//                        continue;
//                    }
//                    else
//                    {
//                        x1 = (-b + sqrt(delta)) / (2 * a);
//                        x2 = (-b - sqrt(delta)) / (2 * a);
//                    }

//                    if (k > 0)
//                    {

//                        if (lwPolyLineMsg->ctrlVector[i].crown > 0)
//                        {
//                            org.X = (x1 > x2) ? x2 : x1;
//                        }
//                        else
//                        {
//                            org.X = (x1 > x2) ? x1 : x2;
//                        }
//                        org.Y = _k * org.X + _b;
//                    }
//                    else
//                    {
//                        if (lwPolyLineMsg->ctrlVector[i].crown < 0)
//                        {
//                            org.X = (x1 > x2) ? x2 : x1;
//                        }
//                        else
//                        {
//                            org.X = (x1 > x2) ? x1 : x2;
//                        }
//                        org.Y = _k * org.X + _b;
//                    }
//                }
//            }

//            /* 求解圆弧夹角A */
//            if ((org.X) == (lwPolyLineMsg->ctrlVector[i].X))
//            {
//                if ((org.Y) > (lwPolyLineMsg->ctrlVector[i].Y))
//                {
//                    angleA = 270;
//                }
//                else
//                {
//                    angleA = 90;
//                }
//            }
//            else if ((org.Y) == (lwPolyLineMsg->ctrlVector[i].Y))
//            {
//                if ((org.X) > (lwPolyLineMsg->ctrlVector[i].X))
//                {
//                    angleA = 180;
//                }
//                else
//                {
//                    angleA = 0;
//                }
//            }
//            else
//            {
//                angleA = fabs(180 * asin(fabs(org.Y - lwPolyLineMsg->ctrlVector[i].Y) / Radius) / M_PI);
//                if ((org.X) < (lwPolyLineMsg->ctrlVector[i].X) && (org.Y) > (lwPolyLineMsg->ctrlVector[i].Y))
//                {
//                    angleA = 360 - angleA;
//                }
//                else if ((org.X) > (lwPolyLineMsg->ctrlVector[i].X) && (org.Y) < (lwPolyLineMsg->ctrlVector[i].Y))
//                {
//                    angleA = 180 - angleA;
//                }
//                else if ((org.X) < (lwPolyLineMsg->ctrlVector[i].X) && (org.Y) < (lwPolyLineMsg->ctrlVector[i].Y))
//                {
//                    //                    angleA = angleA;
//                }
//                else
//                {
//                    angleA = 180 + angleA;
//                }
//            }

//            /* 求解圆弧夹角B */
//            if ((org.X) == (PointNext.X))
//            {
//                if ((org.Y) > (PointNext.Y))
//                {
//                    angleB = 270;
//                }
//                else
//                {
//                    angleB = 90;
//                }
//            }
//            else if ((org.Y) == (PointNext.Y))
//            {
//                if ((org.X) > (PointNext.X))
//                {
//                    angleB = 180;
//                }
//                else
//                {
//                    angleB = 0;
//                }
//            }
//            else
//            {
//                angleB = fabs(180 * asin(fabs(org.Y - PointNext.Y) / Radius) / M_PI);
//                if ((org.X) < (PointNext.X) && (org.Y) > (PointNext.Y))
//                {
//                    angleB = 360 - angleB;
//                }
//                else if ((org.X) > (PointNext.X) && (org.Y) < (PointNext.Y))
//                {
//                    angleB = 180 - angleB;
//                }
//                else if ((org.X) < (PointNext.X) && (org.Y) < (PointNext.Y))
//                {
//                    //                    angleB = angleB;
//                }
//                else
//                {
//                    angleB = 180 + angleB;
//                }
//            }
//            /* 绘制圆弧 */
//            if (lwPolyLineMsg->ctrlVector[i].crown > 0)
//            {
////                msgTmp = GetMinArea_Arc(org.X, org.Y, Radius, angleA, angleB);
//            }
//            else
//            {
////                msgTmp = GetMinArea_Arc(org.X, org.Y, Radius, angleB, angleA);
//            }
//        }

//        if (msg.leftPos.X > msgTmp.leftPos.X)
//        {
//            msg.leftPos = msgTmp.leftPos;
//        }
//        if (msg.rightPos.X < msgTmp.rightPos.X)
//        {
//            msg.rightPos = msgTmp.rightPos;
//        }
//        if (msg.topPos.Y < msgTmp.topPos.Y)
//        {
//            msg.topPos = msgTmp.topPos;
//        }
//        if (msg.bottomPos.Y > msgTmp.bottomPos.Y)
//        {
//            msg.bottomPos = msgTmp.bottomPos;
//        }
//        msg.minRect.X = msg.leftPos.X;
//        msg.minRect.Y = msg.topPos.Y;
//        msg.minRect.Width = msg.rightPos.X - msg.leftPos.X;
//        msg.minRect.Height = msg.topPos.Y - msg.bottomPos.Y;

//        msg.midPos.X = (msg.leftPos.X + msg.rightPos.X) / 2;
//        msg.midPos.Y = (msg.topPos.Y + msg.bottomPos.Y) / 2;
//    }
//    return msg;
//}


/**
 * @brief 求贝塞尔曲线的最小加速度所在t时刻的信息
 * @param {REAL} tStep
 * @param {uint8_t} level
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @return {*}
 * @note
 */
SPLINE_VALUE_t MainWindow::Draw_BezierLine_GetMinAccel(uint8_t level, REAL tStep,
                                                       COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum)
{
    SPLINE_VALUE_t msg = {0, 0};
    COORDINATE_t p00 = {0, 0}, p10 = {0, 0}, p20 = {0, 0};
    REAL accel = 0, accelOld = 0;

    if (ctrlPointNum != level + 1 || level > 3)
    {
        return msg;
    }
    switch (level)
    {
    case 3:
        p00.X = ctrlPointVector[1].X - ctrlPointVector[0].X;
        p00.Y = ctrlPointVector[1].Y - ctrlPointVector[0].Y;

        p10.X = ctrlPointVector[2].X - ctrlPointVector[1].X;
        p10.Y = ctrlPointVector[2].Y - ctrlPointVector[1].Y;

        p20.X = ctrlPointVector[3].X - ctrlPointVector[2].X;
        p20.Y = ctrlPointVector[3].Y - ctrlPointVector[2].Y;

        p00.X = (level - 1) * (p10.X - p00.X);
        p00.Y = (level - 1) * (p10.Y - p00.Y);

        p10.X = (level - 1) * (p20.X - p10.X);
        p10.Y = (level - 1) * (p20.Y - p10.Y);

        accel = sqrt(p00.X * p00.X + p00.Y * p00.Y);
        accelOld = accel;
        /* de Casteljau递推算法 */
        for (REAL t = 0; t < 1; t += tStep)
        {
            COORDINATE_t pxx = {0, 0};
            pxx.X = (level - 2) * ((1 - t) * p00.X + t * p10.X);
            pxx.Y = (level - 2) * ((1 - t) * p00.Y + t * p10.Y);
            accel = sqrt(pxx.X * pxx.X + pxx.Y * pxx.Y);

            if (accel < accelOld)
            {
                accelOld = accel;
                msg.t = t;
                msg.value = accelOld;
            }
        }
        break;
    default:
        break;
    }
    return msg;
}

/**
 * @brief 绘制贝塞尔曲线的一阶导数图像
 * @param {REAL} tStep
 * @param {uint8_t} level
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @return {*}
 * @note
 */
void MainWindow::Draw_BezierLine_FirstDerivative(uint8_t level, REAL tStep,
                                                 COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum)
{
    COORDINATE_t p00 = {0, 0}, p10 = {0, 0}, p20 = {0, 0};

    QPainter painter(&image); //选入绘图设备中。

    if (ctrlPointNum != level + 1 || level > 3)
    {
        return;
    }
    painter.setPen(Qt::blue);
    switch (level)
    {
    case 3:
        p00.X = ctrlPointVector[1].X - ctrlPointVector[0].X;
        p00.Y = ctrlPointVector[1].Y - ctrlPointVector[0].Y;

        p10.X = ctrlPointVector[2].X - ctrlPointVector[1].X;
        p10.Y = ctrlPointVector[2].Y - ctrlPointVector[1].Y;

        p20.X = ctrlPointVector[3].X - ctrlPointVector[2].X;
        p20.Y = ctrlPointVector[3].Y - ctrlPointVector[2].Y;
        /* de Casteljau递推算法 */
        for (REAL t = 0; t < 1; t += tStep)
        {
            COORDINATE_t pxx = {0, 0};
            pxx.X = (level - 1) * ((1 - t) * (1 - t) * p00.X + 2 * t * (1 - t) * p10.X + t * t * p20.X);
            pxx.Y = (level - 1) * ((1 - t) * (1 - t) * p00.Y + 2 * t * (1 - t) * p10.Y + t * t * p20.Y);
            //        LCD_PutPixel((int16_t)p03.X, (int16_t)p03.Y);
            painter.drawPoint(QPoint(pxx.X, pxx.Y));
        }
        break;
    case 2:
        p00.X = ctrlPointVector[1].X - ctrlPointVector[0].X;
        p00.Y = ctrlPointVector[1].Y - ctrlPointVector[0].Y;

        p10.X = ctrlPointVector[2].X - ctrlPointVector[1].X;
        p10.Y = ctrlPointVector[2].Y - ctrlPointVector[1].Y;
        /* de Casteljau递推算法 */
        for (REAL t = 0; t < 1; t += tStep)
        {
            COORDINATE_t pxx = {0, 0};
            pxx.X = (level - 1) * ((1 - t) * p00.X + t * p10.X);
            pxx.Y = (level - 1) * ((1 - t) * p00.Y + t * p10.Y);
            //        LCD_PutPixel((int16_t)p03.X, (int16_t)p03.Y);
            painter.drawPoint(QPoint(pxx.X, pxx.Y));
        }
        break;
    default:
        break;
    }

    painter.end();
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
}

/**
 * @brief 求贝塞尔曲线的最大速度所在t时刻的信息
 * @param {REAL} tStep
 * @param {uint8_t} level
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @return {*}
 * @note
 */
SPLINE_VALUE_t MainWindow::Draw_BezierLine_GetMaxSpeed(uint8_t level, REAL tStep,
                                                       COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum)
{
    SPLINE_VALUE_t msg = {0, 0};
    COORDINATE_t p00 = {0, 0}, p10 = {0, 0}, p20 = {0, 0};
    COORDINATE_t pxx = {0, 0};
    if (ctrlPointNum != level + 1 || level > 3)
    {
        return msg;
    }
    switch (level)
    {
    case 3:
        msg = Draw_BezierLine_GetMinAccel(level, tStep, ctrlPointVector, ctrlPointNum);
        p00.X = ctrlPointVector[1].X - ctrlPointVector[0].X;
        p00.Y = ctrlPointVector[1].Y - ctrlPointVector[0].Y;

        p10.X = ctrlPointVector[2].X - ctrlPointVector[1].X;
        p10.Y = ctrlPointVector[2].Y - ctrlPointVector[1].Y;

        p20.X = ctrlPointVector[3].X - ctrlPointVector[2].X;
        p20.Y = ctrlPointVector[3].Y - ctrlPointVector[2].Y;

        /* de Casteljau递推算法 */
        pxx.X = (level - 1) * ((1 - msg.t) * (1 - msg.t) * p00.X + 2 * msg.t * (1 - msg.t) * p10.X + msg.t * msg.t * p20.X);
        pxx.Y = (level - 1) * ((1 - msg.t) * (1 - msg.t) * p00.Y + 2 * msg.t * (1 - msg.t) * p10.Y + msg.t * msg.t * p20.Y);
        msg.value = sqrt(pxx.X * pxx.X + pxx.Y * pxx.Y);
        break;
    case 2:
        REAL speedTmp;
        p00.X = ctrlPointVector[1].X - ctrlPointVector[0].X;
        p00.Y = ctrlPointVector[1].Y - ctrlPointVector[0].Y;

        p10.X = ctrlPointVector[2].X - ctrlPointVector[1].X;
        p10.Y = ctrlPointVector[2].Y - ctrlPointVector[1].Y;
        msg.value = sqrt(p00.X * p00.X + p00.Y * p00.Y);
        speedTmp = sqrt(p10.X * p10.X + p10.Y * p10.Y);
        if (msg.value < speedTmp)
        {
            msg.value = speedTmp;
            msg.t = 1;
        }
        break;
    default:
        break;
    }

    return msg;
}

/**
 * @brief 绘制贝塞尔曲线-最大支持3阶
 * @param {uint8_t} level
 * @param {REAL} tStep t步长
 * @param {REAL} ratio 坐标最小精度
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @return {*}
 * @note
 */
MIN_AREA_T MainWindow::Draw_BezierLine(uint8_t level, REAL tStep, REAL ratio,
                                       COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum)
{
    uint16_t maxLen = 0; //用于机选自适应步长的中间值
    COORDINATE_t p00 = {0, 0}, p10 = {0, 0}, p20 = {0, 0}, p30 = {0, 0};
    QPainter painter(&image); //选入绘图设备中。

    MIN_AREA_T msg = {
        .minRect = {ctrlPointVector[0].X, ctrlPointVector[0].Y, 0, 0},
        .leftPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .rightPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .topPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .bottomPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y}};

    if (ctrlPointNum != level + 1 || level > 3)
    {
        memset(&msg, 0, sizeof(MIN_AREA_T));
        return msg;
    }

    /* 计算分辨率 */
    if (tStep > 1 || tStep <= 0)
    {
        SPLINE_VALUE_t accelMsg = {0, 0};
        tStep = 0;
        /* 获取加速度最小，速度最大时刻的t值与速度值 */
        accelMsg = Draw_BezierLine_GetMaxSpeed(level, ratio / 100, ctrlPointVector, ctrlPointNum);

        /* 获取最小分辨率的数值的向上取整 */
        maxLen = (int)(accelMsg.value * 1.414 * 2 / ratio);
        while (pow(10, tStep) < maxLen)
        {
            tStep++;
        }
        /* 获取最小步长 */
        tStep = 1 / (pow(10, tStep));
    }
    /* 显示调试用控制点 */
    if (BEZIER_SHOW_CTRL)
    {
        painter.setPen(Qt::red);
        for (uint8_t i = 0; i < ctrlPointNum - 1; i++)
        {
            QLineF line1(ctrlPointVector[i].X, ctrlPointVector[i].Y,
                         ctrlPointVector[i + 1].X, ctrlPointVector[i + 1].Y);
            painter.drawLine(line1);
        }
    }

    painter.setPen(Qt::blue);
    switch (level)
    {
    case 3:
        p00 = ctrlPointVector[0];
        p10 = ctrlPointVector[1];
        p20 = ctrlPointVector[2];
        p30 = ctrlPointVector[3];
        /* de Casteljau递推算法 */
        for (REAL t = 0; t < 1; t += tStep)
        {
            COORDINATE_t p01 = {0, 0}, p11 = {0, 0}, p21 = {0, 0}; //一阶节点
            COORDINATE_t p02 = {0, 0}, p12 = {0, 0};               //二阶节点
            COORDINATE_t p03 = {0, 0};                             //三阶节点
            /* 一次递推 */
            p01.X = (1 - t) * p00.X + t * p10.X;
            p11.X = (1 - t) * p10.X + t * p20.X;
            p21.X = (1 - t) * p20.X + t * p30.X;
            p01.Y = (1 - t) * p00.Y + t * p10.Y;
            p11.Y = (1 - t) * p10.Y + t * p20.Y;
            p21.Y = (1 - t) * p20.Y + t * p30.Y;
            /* 二次递推 */
            p02.X = (1 - t) * p01.X + t * p11.X;
            p12.X = (1 - t) * p11.X + t * p21.X;
            p02.Y = (1 - t) * p01.Y + t * p11.Y;
            p12.Y = (1 - t) * p11.Y + t * p21.Y;
            /* 三次递推 */
            p03.X = (1 - t) * p02.X + t * p12.X;
            p03.Y = (1 - t) * p02.Y + t * p12.Y;
            /* 获取最左边点 */
            if (p03.X < msg.leftPos.X)
            {
                msg.leftPos.X = p03.X;
                msg.leftPos.Y = p03.Y;
            }
            /* 获取最右边点 */
            if (p03.X > msg.rightPos.X)
            {
                msg.rightPos.X = p03.X;
                msg.rightPos.Y = p03.Y;
            }
            /* 获取最上边点 */
            if (p03.Y > msg.topPos.Y)
            {
                msg.topPos.X = p03.X;
                msg.topPos.Y = p03.Y;
            }
            /* 获取最下边点 */
            if (p03.Y < msg.bottomPos.Y)
            {
                msg.bottomPos.X = p03.X;
                msg.bottomPos.Y = p03.Y;
            }

            //        LCD_PutPixel((int16_t)p03.X, (int16_t)p03.Y);
            painter.drawPoint(QPoint(p03.X, p03.Y));
        }
        break;
    case 2:
        p00 = ctrlPointVector[0];
        p10 = ctrlPointVector[1];
        p20 = ctrlPointVector[2];
        /* de Casteljau递推算法 */
        for (REAL t = 0; t < 1; t += tStep)
        {
            COORDINATE_t p01 = {0, 0}, p11 = {0, 0}; //一阶节点
            COORDINATE_t p02 = {0, 0};               //二阶节点
            /* 一次递推 */
            p01.X = (1 - t) * p00.X + t * p10.X;
            p11.X = (1 - t) * p10.X + t * p20.X;
            p01.Y = (1 - t) * p00.Y + t * p10.Y;
            p11.Y = (1 - t) * p10.Y + t * p20.Y;
            /* 二次递推 */
            p02.X = (1 - t) * p01.X + t * p11.X;
            p02.Y = (1 - t) * p01.Y + t * p11.Y;

            /* 获取最左边点 */
            if (p02.X < msg.leftPos.X)
            {
                msg.leftPos.X = p02.X;
                msg.leftPos.Y = p02.Y;
            }
            /* 获取最右边点 */
            if (p02.X > msg.rightPos.X)
            {
                msg.rightPos.X = p02.X;
                msg.rightPos.Y = p02.Y;
            }
            /* 获取最上边点 */
            if (p02.Y > msg.topPos.Y)
            {
                msg.topPos.X = p02.X;
                msg.topPos.Y = p02.Y;
            }
            /* 获取最下边点 */
            if (p02.Y < msg.bottomPos.Y)
            {
                msg.bottomPos.X = p02.X;
                msg.bottomPos.Y = p02.Y;
            }

            //        LCD_PutPixel((int16_t)p03.X, (int16_t)p03.Y);
            painter.drawPoint(QPoint(p02.X, p02.Y));
        }
        break;
    case 1:
        break;
    default:
        break;
    }

    msg.minRect.X = msg.leftPos.X;
    msg.minRect.Y = msg.bottomPos.Y;
    msg.minRect.Width = msg.rightPos.X - msg.leftPos.X;
    msg.minRect.Height = msg.topPos.Y - msg.bottomPos.Y;

    /* 显示最小矩形 */
    if (BEZIER_SHOW_CTRL)
    {
    }
    painter.end();
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
    return msg;
}

/**
 * @brief B—spline 的基函数系数递推公式,逆推法,递归
 * @param {uint16_t} level
 * @param {REAL} *KnotVector
 * @param {uint16_t} KnotNum
 * @param {uint16_t} KnotIndex
 * @param {REAL} t
 * @return {*}
 * @note
 */
REAL MainWindow::Draw_Bspline_GetBidegN(uint16_t level, REAL *KnotVector, uint16_t KnotNum,
                                        uint16_t KnotIndex, REAL t)
{
    REAL Bideg = 0;    //基本函数值
    REAL Knot_Max = 1; //最大节点值
    /* 获取最大的节点值 */
    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (Knot_Max < KnotVector[i])
        {
            Knot_Max = KnotVector[i];
        }
    }

    /* 标准化处理 */
    // for (uint16_t i = 0; i <= KnotNum; ++i)
    // {
    // 	KnotVector[i] /= Knot_Max;
    // }
    // Knot_Max = 1;

    if (level != 0)
    {
        REAL aA_a = Knot_Max * t - KnotVector[KnotIndex];                          //第一项分子
        REAL aA_A = KnotVector[KnotIndex + level] - KnotVector[KnotIndex];         //第一项分母
        REAL bB_b = KnotVector[KnotIndex + level + 1] - Knot_Max * t;              //第二项分子
        REAL bB_B = KnotVector[KnotIndex + level + 1] - KnotVector[KnotIndex + 1]; //第二项分母
        REAL aA_aA;
        REAL bB_bB;
        /* 规定0/0=0 */
        if (aA_a == 0 || aA_A == 0)
        {
            aA_aA = 0;
        }
        else
        {
            aA_aA = aA_a / aA_A;
        }
        if (bB_b == 0 || bB_B == 0)
        {
            bB_bB = 0;
        }
        else
        {
            bB_bB = bB_b / bB_B;
        }
        Bideg = aA_aA * Draw_Bspline_GetBidegN(level - 1, KnotVector, KnotNum, KnotIndex, t) +
                bB_bB * Draw_Bspline_GetBidegN(level - 1, KnotVector, KnotNum, KnotIndex + 1, t);
    }
    else
    {
        if (Knot_Max * t >= KnotVector[KnotIndex] && Knot_Max * t <= KnotVector[KnotIndex + 1])
        {
            Bideg = 1;
        }
        else
        {
            Bideg = 0;
        }
    }

    return Bideg;
}

/**
 * @brief B—spline 的基函数系数递推公式,正推法
 * @param {uint16_t} level
 * @param {REAL} *KnotVector
 * @param {uint16_t} KnotNum
 * @param {uint16_t} KnotIndex
 * @param {REAL} t
 * @return {*}
 * @note
 */
REAL MainWindow::Draw_Bspline_GetBidegP(uint16_t level, REAL *KnotVector, uint16_t KnotNum,
                                        uint16_t KnotIndex, REAL t)
{
    REAL Bideg = 0;                  //基本函数值
    REAL Bi_n[KNOT_NUM_LIMIT] = {0}; //基本函数值
    REAL Knot_Max = 1;               //最大节点值
    /* 获取最大的节点值 */
    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (Knot_Max < KnotVector[i])
        {
            Knot_Max = KnotVector[i];
        }
    }

    /* 标准化处理 */
    // for (uint16_t i = 0; i <= KnotNum; ++i)
    // {
    // 	KnotVector[i] /= Knot_Max;
    // }
    // Knot_Max = 1;

    for (uint16_t levelTmp = 0; levelTmp <= level; levelTmp++)
    {
        if (levelTmp != 0)
        {
            REAL aA_a = Knot_Max * t - KnotVector[KnotIndex];                             //第一项分子
            REAL aA_A = KnotVector[KnotIndex + levelTmp] - KnotVector[KnotIndex];         //第一项分母
            REAL bB_b = KnotVector[KnotIndex + levelTmp + 1] - Knot_Max * t;              //第二项分子
            REAL bB_B = KnotVector[KnotIndex + levelTmp + 1] - KnotVector[KnotIndex + 1]; //第二项分母
            REAL aA_aA;
            REAL bB_bB;
            /* 规定0/0=0 */
            if (aA_a == 0 || aA_A == 0)
            {
                aA_aA = 0;
            }
            else
            {
                aA_aA = aA_a / aA_A;
            }
            if (bB_b == 0 || bB_B == 0)
            {
                bB_bB = 0;
            }
            else
            {
                bB_bB = bB_b / bB_B;
            }
            Bi_n[levelTmp] = aA_aA * Bi_n[levelTmp - 1] + bB_bB * Bi_n[levelTmp - 1];
            Bideg = Bi_n[levelTmp];
        }
        else
        {

            if (Knot_Max * t >= KnotVector[KnotIndex] && Knot_Max * t <= KnotVector[KnotIndex + 1])
            {
                Bi_n[0] = 1;
            }
            else
            {
                Bi_n[0] = 0;
            }
            if (level == 0)
            {
                Bideg = Bi_n[0];
            }
        }
    }

    return Bideg;
}

/**
 * @brief 获取B-spline的阶段数组-Clamped列表法/均匀参数化法
 * @param {uint16_t} level
 * @param {uint16_t} ctrlPointNum
 * @param {REAL} *KnotVector
 * @return {uint16_t} 节点矢量数量
 * @note
 */
uint16_t MainWindow::Draw_Bspline_CreateKnotVectorClamped(uint16_t level, uint16_t ctrlPointNum, REAL *KnotVector)
{
    /* 计算节点值，使用非均匀Clamped的方法 前阶数+1个数的节点值 = 0  后阶数+1个数的节点值 = 1 其余为0-1的均分 */
    uint16_t KnotNum = ctrlPointNum + level + 1;         //节点数
    uint16_t KnotStep = 1;                               //归一化处理的步长
    uint16_t KnotNeedAveNum = KnotNum - (level + 1) * 2; //需要均分的节点数

    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (i <= level)
        {
            KnotVector[i] = 0;
        }
        else if (i > level && i < KnotNum - level)
        {
            KnotVector[i] = KnotVector[i - 1] + KnotStep;
        }
        else if (i >= KnotNum - level)
        {
            KnotVector[i] = KnotVector[i - 1];
        }
    }
    return KnotNum;
}

/**
 * @brief 获取B-spline的阶段数组-累计弦长法
 * @param {uint16_t} level
 * @param {uint16_t} ctrlPointNum
 * @param {REAL} *KnotVector
 * @return {uint16_t} 节点矢量数量
 * @note
 */
uint16_t MainWindow::Draw_Bspline_CreateKnotVectorSumChord(uint16_t level,
                                                           COORDINATE_t *fitPointVector, uint16_t fitPointNum, REAL *KnotVector)
{
    /* 计算节点值，使用非均匀Clamped的方法 前阶数+1个数的节点值 = 0  后阶数+1个数的节点值 = 1 其余为0-1的均分 */
    uint16_t KnotNum = fitPointNum + 2 + level + 1; //节点数
    uint16_t KnotStep = 1;                          //归一化处理的步长

    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (i <= level)
        {
            KnotVector[i] = 0;
        }

        else if (i > level && i < KnotNum - level)
        {
            KnotVector[i] = sqrt((fitPointVector[i - level].X - fitPointVector[i - level - 1].X) * (fitPointVector[i - level].X - fitPointVector[i - level - 1].X) + (fitPointVector[i - level].Y - fitPointVector[i - level - 1].Y) * (fitPointVector[i - level].Y - fitPointVector[i - level - 1].Y)) +
                    KnotVector[i - 1];
        }

        else if (i >= KnotNum - level)
        {
            KnotVector[i] = KnotVector[i - 1];
        }
    }
    return KnotNum;
}

COORDINATE_t MainWindow::Draw_Bspline_DeBoor(uint16_t level,  REAL t,
                                             COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum,
                                             REAL *KnotVector, uint16_t KnotNum)
{
    COORDINATE_t Pxx = {0, 0};       //拟合点
    REAL Knot_Max = 0;               //最大节点值
    uint16_t i = 0;
    uint16_t r = level - 1;
    uint16_t l = level - r;
    REAL an_1[3] = {0};
    REAL an_2[2] = {0};
    REAL an_3[1] = {0};
    COORDINATE_t dn_1[3];
    COORDINATE_t dn_2[2];
    COORDINATE_t dn_3[1];


    /* 判断节点数是否达标 节点数 = 控制点数 + 阶数 + 1 */
    if (KnotNum != ctrlPointNum + level + 1)
    {
        return Pxx;
    }
    /* 获取最大的节点值 */
    for (uint16_t j = 0; j < KnotNum; ++j)
    {
        if (Knot_Max < KnotVector[j])
        {
            Knot_Max = KnotVector[j];
        }
    }
    /* de boor递推算法 */
    Pxx.X = 0;
    Pxx.Y = 0;

    t *= Knot_Max;
    for(;i < KnotNum - 1;i++)
    {
        if(t >= KnotVector[i] && t <KnotVector[i + 1])
        {
            break;
        }
    }
    if(level == 3)
    {
        r = level - 1;
        l = level - r;
        for(uint16_t j = i - level;j <= i - level + r;j++)
        {
            if((t - KnotVector[j + l] == 0)||(KnotVector[j + level + 1] - KnotVector[j + l] == 0))
            {
                an_1[j - i +level] = 0;
            }
            else
            {
                an_1[j - i +level] = (t - KnotVector[j + l])/(KnotVector[j + level + 1] - KnotVector[j + l]);
            }
            dn_1[j - i +level].X = (1 - an_1[j - i +level])*ctrlPointVector[j].X + an_1[j - i +level]*ctrlPointVector[j+1].X;
            dn_1[j - i +level].Y = (1 - an_1[j - i +level])*ctrlPointVector[j].Y + an_1[j - i +level]*ctrlPointVector[j+1].Y;
        }

        r = level - 2;
        l = level - r;
        for(uint16_t j = i - level;j <= i - level + r;j++)
        {
            if((t - KnotVector[j + l] == 0)||(KnotVector[j + level + 1] - KnotVector[j + l] == 0))
            {
                an_2[j - i +level] = 0;
            }
            else
            {
                an_2[j - i +level] = (t - KnotVector[j + l])/(KnotVector[j + level + 1] - KnotVector[j + l]);
            }
            dn_2[j - i +level].X = (1 - an_2[j - i +level])*dn_1[j - i +level].X + an_2[j - i +level]*dn_1[j+1 - i +level].X;
            dn_2[j - i +level].Y = (1 - an_2[j - i +level])*dn_1[j - i +level].Y + an_2[j - i +level]*dn_1[j+1 - i +level].Y;
        }

        r = level - 3;
        l = level - r;
        for(uint16_t j = i - level;j <= i - level + r;j++)
        {
            if((t - KnotVector[j + l] == 0)||(KnotVector[j + level + 1] - KnotVector[j + l] == 0))
            {
                an_3[j - i +level] = 0;
            }
            else
            {
                an_3[j - i +level] = (t - KnotVector[j + l])/(KnotVector[j + level + 1] - KnotVector[j + l]);
            }
            dn_3[j - i +level].X = (1 - an_3[j - i +level])*dn_2[j - i +level].X + an_3[j - i +level]*dn_2[j+1 - i +level].X;
            dn_3[j - i +level].Y = (1 - an_3[j - i +level])*dn_2[j - i +level].Y + an_3[j - i +level]*dn_2[j+1 - i +level].Y;
        }
    }
    else if(level == 2)
    {
        r = level - 1;
        l = level - r;
        for(uint16_t j = i - level;j <= i - level + r;j++)
        {
            if((t - KnotVector[j + l] == 0)||(KnotVector[j + level + 1] - KnotVector[j + l] == 0))
            {
                an_2[j - i +level] = 0;
            }
            else
            {
                an_2[j - i +level] = (t - KnotVector[j + l])/(KnotVector[j + level + 1] - KnotVector[j + l]);
            }
            dn_2[j - i +level].X = (1 - an_2[j - i +level])*ctrlPointVector[j].X + an_2[j - i +level]*ctrlPointVector[j+1].X;
            dn_2[j - i +level].Y = (1 - an_2[j - i +level])*ctrlPointVector[j].Y + an_2[j - i +level]*ctrlPointVector[j+1].Y;
        }

        r = level - 2;
        l = level - r;
        for(uint16_t j = i - level;j <= i - level + r;j++)
        {
            if((t - KnotVector[j + l] == 0)||(KnotVector[j + level + 1] - KnotVector[j + l] == 0))
            {
                an_3[j - i +level] = 0;
            }
            else
            {
                an_3[j - i +level] = (t - KnotVector[j + l])/(KnotVector[j + level + 1] - KnotVector[j + l]);
            }
            dn_3[j - i +level].X = (1 - an_3[j - i +level])*dn_2[j - i +level].X + an_3[j - i +level]*dn_2[j+1 - i +level].X;
            dn_3[j - i +level].Y = (1 - an_3[j - i +level])*dn_2[j - i +level].Y + an_3[j - i +level]*dn_2[j+1 - i +level].Y;
        }
    }
    Pxx.X = dn_3[0].X;
    Pxx.Y = dn_3[0].Y;
    return  Pxx;
}

COORDINATE_t MainWindow::Draw_Bspline_DeBoor_FirstDerivative(uint16_t level,  REAL t,
                                                             COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum,
                                                             REAL *KnotVector, uint16_t KnotNum)
{
    COORDINATE_t Pxx = {0, 0};       //拟合点
    REAL Knot_Max = 0;               //最大节点值
    uint16_t i = 0;
    uint16_t r = 1;
    uint16_t l = 1;
    REAL an_1[3] = {0};
    REAL an_2[2] = {0};
    REAL an_3[1] = {0};
    COORDINATE_t dn_1[3];
    COORDINATE_t dn_2[2];
    COORDINATE_t dn_3[1];


    /* 判断节点数是否达标 节点数 = 控制点数 + 阶数 + 1 */
    if (KnotNum != ctrlPointNum + level + 1)
    {
        return Pxx;
    }
    if(level > 3)
    {
        return Pxx;
    }
    /* 获取最大的节点值 */
    for (uint16_t j = 0; j < KnotNum; ++j)
    {
        if (Knot_Max < KnotVector[j])
        {
            Knot_Max = KnotVector[j];
        }
    }

    //    for (uint16_t j = 0; j < KnotNum; ++j)
    //    {
    //        KnotVector[j] /= Knot_Max;
    //    }
    /* de boor递推算法 */
    Pxx.X = 0;
    Pxx.Y = 0;

    t *= Knot_Max;
    for(;i < KnotNum - 1;i++)
    {
        if(t >= KnotVector[i] && t <KnotVector[i + 1])
        {
            break;
        }
    }

    for(l = 1;l <= r; l++)
    {
        for(uint16_t j = i-level;j<= i-l;j++)
        {
            REAL beta;
            if(KnotVector[j + level + 1] - KnotVector[j + 1] == 0)
            {
                beta = 0;
            }
            else
            {
                beta = (level - l + 1)/(KnotVector[j + level + 1] - KnotVector[j + 1]);
                beta *= Knot_Max;
            }
            dn_1[j - i + level].X = beta * (ctrlPointVector[j - i + level + 1].X - ctrlPointVector[j - i + level].X);
            dn_1[j - i + level].Y = beta * (ctrlPointVector[j - i + level + 1].Y - ctrlPointVector[j - i + level].Y);
        }
    }

    for(l = 1;l <= level - r; l++)
    {
        for(uint16_t j = i-level;j<= i-l-r;j++)
        {
            REAL alpha,du;
            du = KnotVector[j + level + 1] - KnotVector[j + r + 1];
            if(du == 0)
            {
                alpha = 0;
            }
            else
            {
                alpha = (t - KnotVector[j + r + 1])/(du);
            }
            dn_3[j - i + level].X = (1-alpha) * dn_1[j - i + level].X + alpha*dn_1[j - i + level + 1].X;
            dn_3[j - i + level].Y = (1-alpha) * dn_1[j - i + level].Y + alpha*dn_1[j - i + level + 1].Y;
        }
    }

    Pxx.X = dn_3[0].X;
    Pxx.Y = dn_3[0].Y;
    return  Pxx;
}


/**
 * @brief 获取B样条曲线最大速度所在t时刻的信息-需要输入节点矢量表
 * @param {REAL} tStep
 * @param {uint8_t} level
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @param {REAL} *KnotVector
 * @param {uint16_t} KnotNum
 * @return {*}
 * @note
 */
SPLINE_VALUE_t MainWindow::Draw_Bspline_GetMaxSpeedWithKnot(uint8_t level, REAL tStep,
                                                            COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum,
                                                            REAL *KnotVector, uint16_t KnotNum)
{

    COORDINATE_t Pxx = {0, 0}; //拟合点
    REAL Knot_Max = 0;         //最大节点值
    REAL maxLen = 0;           //用于机选自适应步长的中间值
    SPLINE_VALUE_t msg = {0, 0};
    REAL speed = 0;


    if (KnotNum != ((ctrlPointNum) + (level) + 1))
    {
        return msg;
    }
    if(level >= 3)
    {


        /* 获取最大的节点值 */
        for (uint16_t i = 0; i < KnotNum; ++i)
        {
            if (Knot_Max < KnotVector[i])
            {
                Knot_Max = KnotVector[i];
            }
        }

        /* de boor递推算法 */
        for (REAL t = 0; t < 1; t += tStep)
        {
            Pxx.X = 0;
            Pxx.Y = 0;
            //        for (uint16_t n = 0; n < ctrlPointNum - 1; n++)
            //        {
            //            REAL Bideg = Draw_Bspline_GetBidegN(level - 1, KnotVector, KnotNum, n + 1, t);
            //            REAL delta = (ctrlPointVector[n + 1].X - ctrlPointVector[n].X);
            //            REAL coefficient = (Knot_Max * level) / (KnotVector[n + level + 1] - KnotVector[n + 1]);
            //            Pxx.X += Bideg * delta * coefficient;
            //            delta = (ctrlPointVector[n + 1].Y - ctrlPointVector[n].Y);
            //            Pxx.Y += Bideg * delta * coefficient;
            //        }
            Pxx = Draw_Bspline_DeBoor_FirstDerivative(level,t,ctrlPointVector,ctrlPointNum,KnotVector,KnotNum);
            speed = sqrt(Pxx.X * Pxx.X + Pxx.Y * Pxx.Y);
            if (speed > msg.value)
            {
                msg.value = speed;
                msg.t = t;
            }
        }
    }
    else
    {
        maxLen = 0;
        REAL maxLenTmp = 0;
        for(uint16_t i = 0; i < ctrlPointNum - 1; i++)
        {
            maxLenTmp = sqrt((ctrlPointVector[i].X - ctrlPointVector[i+1].X) *(ctrlPointVector[i].X - ctrlPointVector[i+1].X)+
                    (ctrlPointVector[i].Y - ctrlPointVector[i+1].Y) *(ctrlPointVector[i].Y - ctrlPointVector[i+1].Y));
            if(maxLenTmp > maxLen)
            {
                maxLen = maxLenTmp;
                msg.t = i/ctrlPointNum;
                msg.value = maxLen;
            }
        }
    }
    return msg;
}

/**
 * @brief 绘制B样条曲线-含节点数组
 * @param {REAL} ratio 坐标最小分辨率
 * @param {uint16_t} level
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @param {REAL} *KnotVector
 * @param {uint16_t} KnotNum
 * @return {*}
 * @note
 */
MIN_AREA_T MainWindow::Draw_BsplineWithKnot(uint16_t level, REAL tStep, REAL ratio,
                                            COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum,
                                            REAL *KnotVector, uint16_t KnotNum)
{

    COORDINATE_t Pxx = {0, 0};       //拟合点
    COORDINATE_t PxxBefore = {0, 0}; //拟合点
    REAL Knot_Max = 0;               //最大节点值
    QPainter painter(&image);        //选入绘图设备中。
    painter.setPen(Qt::blue);
    MIN_AREA_T msg = {
        .minRect = {ctrlPointVector[0].X, ctrlPointVector[0].Y, 0, 0},
        .leftPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .rightPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .topPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .bottomPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y}};
    /* 判断节点数是否达标 节点数 = 控制点数 + 阶数 + 1 */
    if (KnotNum != ctrlPointNum + level + 1)
    {
        return msg;
    }
    /* 获取最大的节点值 */
    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (Knot_Max < KnotVector[i])
        {
            Knot_Max = KnotVector[i];
        }
    }

    /* 计算分辨率 */
    if (tStep > 1 || tStep <= 0)
    {
        REAL maxLen = 0; //用于机选自适应步长的中间值
        SPLINE_VALUE_t splineMsg = {0, 0};
        tStep = 1;
        splineMsg = Draw_Bspline_GetMaxSpeedWithKnot(level, ratio / 100, ctrlPointVector, ctrlPointNum, KnotVector, KnotNum);
        /* 获取最小分辨率的数值的向上取整 */
        maxLen = (int)(splineMsg.value * 1.414 * 2 / ratio);
        while (pow(10, tStep) < maxLen)
        {
            tStep++;
        }
        /* 获取最小步长 */
        tStep = 1 / (pow(10, tStep));
    }

    /* 绘制曲线 */
    for (REAL t = 0; t < 1; t += tStep)
    {
        Pxx.X = 0;
        Pxx.Y = 0;
        //        for (uint16_t n = 0; n < ctrlPointNum; n++)
        //        {
        //            Pxx.X += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
        //            Pxx.Y += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
        //        }
        Pxx = Draw_Bspline_DeBoor(level,t,ctrlPointVector,ctrlPointNum,KnotVector,KnotNum);
        if (Pxx.X == 0 && Pxx.Y == 0)
        {
            for (uint16_t n = 0; n < ctrlPointNum; n++)
            {
                Pxx.X += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
                Pxx.Y += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
            }
        }
        /* 获取最左边点 */
        if (Pxx.X < msg.leftPos.X)
        {
            msg.leftPos.X = Pxx.X;
            msg.leftPos.Y = Pxx.Y;
        }
        /* 获取最右边点 */
        if (Pxx.X > msg.rightPos.X)
        {
            msg.rightPos.X = Pxx.X;
            msg.rightPos.Y = Pxx.Y;
        }
        /* 获取最上边点 */
        if (Pxx.Y > msg.topPos.Y)
        {
            msg.topPos.X = Pxx.X;
            msg.topPos.Y = Pxx.Y;
        }
        /* 获取最下边点 */
        if (Pxx.Y < msg.bottomPos.Y)
        {
            msg.bottomPos.X = Pxx.X;
            msg.bottomPos.Y = Pxx.Y;
        }

        //        LCD_PutPixel((int16_t)Pxx.X, (int16_t)Pxx.Y);
        PxxBefore.X = Pxx.X;
        PxxBefore.Y = Pxx.Y;
        painter.drawPoint(QPoint(Pxx.X, Pxx.Y));
    }

    msg.minRect.X = msg.leftPos.X;
    msg.minRect.Y = msg.bottomPos.Y;
    msg.minRect.Width = msg.rightPos.X - msg.leftPos.X;
    msg.minRect.Height = msg.topPos.Y - msg.bottomPos.Y;

    /* 显示最小矩形 */
    if (B_SPLINE_SHOW_CTRL)
    {
    }
    painter.end();
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
    return msg;
}

/**
 * @brief 绘制B-spline,不含节点数组,自生产均匀节点数组
 * @param {uint16_t} level
 * @param {REAL} tStep
 * @param {REAL} ratio
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @return {*}
 * @note
 */
MIN_AREA_T MainWindow::Draw_BsplineClamped(uint16_t level, REAL tStep, REAL ratio,
                                           COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum)
{
    REAL KnotVector[KNOT_NUM_LIMIT] = {0}; //节点数组
    COORDINATE_t Pxx = {0, 0};             //拟合点
    REAL Knot_Max = 0;                     //最大节点值
    QPainter painter(&image);              //选入绘图设备中。
    painter.setPen(Qt::blue);
    uint16_t KnotNum = ctrlPointNum + level + 1; //节点数量

    MIN_AREA_T msg = {
        .minRect = {ctrlPointVector[0].X, ctrlPointVector[0].Y, 0, 0},
        .leftPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .rightPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .topPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .bottomPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y}};
    Draw_Bspline_CreateKnotVectorClamped(level, ctrlPointNum, KnotVector);
    /* 获取最大的节点值 */
    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (Knot_Max < KnotVector[i])
        {
            Knot_Max = KnotVector[i];
        }
    }

    /* 计算分辨率 */
    if (tStep > 1 || tStep <= 0)
    {
        SPLINE_VALUE_t splineMsg = {0, 0};
        REAL maxLen = 0; //用于机选自适应步长的中间值
        tStep = 1;
        splineMsg = Draw_Bspline_GetMaxSpeedWithKnot(level, ratio / 100, ctrlPointVector, ctrlPointNum, KnotVector, KnotNum);
        /* 获取最小分辨率的数值的向上取整 */
        maxLen = (int)(splineMsg.value * 1.414 * 2 / ratio);
        while (pow(10, tStep) < maxLen)
        {
            tStep++;
        }
        /* 获取最小步长 */
        tStep = 1 / pow(10, tStep);
    }
    /* de boor递推算法 */
    for (REAL t = 0; t < 1; t += tStep)
    {
        Pxx.X = 0;
        Pxx.Y = 0;
        //        for (uint16_t n = 0; n < ctrlPointNum; n++)
        //        {
        //            Pxx.X += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
        //            Pxx.Y += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
        //        }
        Pxx = Draw_Bspline_DeBoor(level,t,ctrlPointVector,ctrlPointNum,KnotVector,KnotNum);
        /* 获取最左边点 */
        if (Pxx.X < msg.leftPos.X)
        {
            msg.leftPos.X = Pxx.X;
            msg.leftPos.Y = Pxx.Y;
        }
        /* 获取最右边点 */
        if (Pxx.X > msg.rightPos.X)
        {
            msg.rightPos.X = Pxx.X;
            msg.rightPos.Y = Pxx.Y;
        }
        /* 获取最上边点 */
        if (Pxx.Y > msg.topPos.Y)
        {
            msg.topPos.X = Pxx.X;
            msg.topPos.Y = Pxx.Y;
        }
        /* 获取最下边点 */
        if (Pxx.Y < msg.bottomPos.Y)
        {
            msg.bottomPos.X = Pxx.X;
            msg.bottomPos.Y = Pxx.Y;
        }

        //        LCD_PutPixel((int16_t)Pxx.X, (int16_t)Pxx.Y);
        painter.drawPoint(QPoint(Pxx.X, Pxx.Y));
    }

    msg.minRect.X = msg.leftPos.X;
    msg.minRect.Y = msg.bottomPos.Y;
    msg.minRect.Width = msg.rightPos.X - msg.leftPos.X;
    msg.minRect.Height = msg.topPos.Y - msg.bottomPos.Y;

    /* 显示最小矩形 */
    if (B_SPLINE_SHOW_CTRL)
    {
    }
    painter.end();
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
    return msg;
}

/**
 * @brief 绘制B-spline,通过拟合点进行绘制
 * @param {uint16_t} level
 * @param {REAL} tStep
 * @param {REAL} ratio
 * @param {COORDINATE_t} *FitVector
 * @param {uint16_t} FitNum
 * @return {*}
 * @note
 */
MIN_AREA_T MainWindow::Draw_BsplineWithFit(uint16_t level, REAL ratio,
                                           COORDINATE_t *FitVector, uint16_t FitNum)
{
    REAL KnotVector[KNOT_NUM_LIMIT] = {0}; //节点数组
    REAL dKnotVector[KNOT_NUM_LIMIT] = {0};
    COORDINATE_t ctrlPointVector[CTRL_POINT_NUM_LIMIT];
    uint16_t ctrlPointNum = FitNum + 2;
    REAL Knot_Max = 0;                           //最大节点值
    uint16_t KnotNum = ctrlPointNum + level + 1; //节点数量
    MIN_AREA_T msg = {
        .minRect = {ctrlPointVector[0].X, ctrlPointVector[0].Y, 0, 0},
        .leftPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .rightPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .topPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y},
        .bottomPos = {ctrlPointVector[0].X, ctrlPointVector[0].Y}};
    /* 获取节点矢量 */
    //  Draw_Bspline_CreateKnotVectorClamped(level, ctrlPointNum, KnotVector);//
    Draw_Bspline_CreateKnotVectorSumChord(level, FitVector, FitNum, KnotVector);
    /* 获取最大的节点值 */
    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (Knot_Max < KnotVector[i])
        {
            Knot_Max = KnotVector[i];
        }
    }
    /* 获取dKnot */
    for (uint16_t i = level; i < FitNum + level - 1; i++)
    {
        dKnotVector[i] = KnotVector[i + 1] - KnotVector[i];
        dKnotVector[i] /= Knot_Max;
    }

    if (1)
    {

        Matrix Ai;  //节点矢量矩阵
        Matrix Eix; //曲线切矢X
        Matrix Eiy;
        Matrix Dix; //控制点坐标X
        Matrix Diy;
        Matrix AiInverse; //节点矢量矩阵的逆
        Ai = Matrix_CreateZero(FitNum, FitNum);
        Eix = Matrix_CreateZero(FitNum, 1);
        Eiy = Matrix_CreateZero(FitNum, 1);
        Matrix_Zero(Ai);
        Matrix_Zero(Eix);
        Matrix_Zero(Eiy);

        /* 设置边界条件 */
        Ai->data[0][0] = 1;
        Ai->data[FitNum - 1][FitNum - 1] = 1;

        //边界条件1 固定方向
        //            Eix->data[0][0] = FitVector[0].X + dKnotVector[level]/level* 0;
        //            Eiy->data[0][0] = FitVector[0].Y - dKnotVector[level]/level* 1;
        //            Eix->data[FitNum-1][0] = FitVector[FitNum-1].X + dKnotVector[FitNum + 1]/level* (-1);
        //            Eiy->data[FitNum-1][0] = FitVector[FitNum-1].Y - dKnotVector[FitNum + 1]/level* 0;

        //边界条件2 e0 = p0 ; en-1 = pn-1
        //        Eix->data[0][0] = FitVector[0].X;
        //        Eiy->data[0][0] = FitVector[0].Y;
        //        Eix->data[FitNum - 1][0] = FitVector[FitNum - 1].X;
        //        Eiy->data[FitNum - 1][0] = FitVector[FitNum - 1].Y;

        //边界条件3
        Eix->data[0][0] = FitVector[0].X + dKnotVector[level]/level* (FitVector[1].X - FitVector[0].X);
        Eiy->data[0][0] = FitVector[0].Y + dKnotVector[level]/level* (FitVector[1].Y - FitVector[0].Y);
        Eix->data[FitNum-1][0] = FitVector[FitNum-1].X - dKnotVector[FitNum + 1]/level* (FitVector[FitNum-1].X - FitVector[FitNum-2].X);
        Eiy->data[FitNum-1][0] = FitVector[FitNum-1].Y - dKnotVector[FitNum + 1]/level* (FitVector[FitNum-1].Y - FitVector[FitNum-2].Y);

        /* 节点矢量表计算 */
        Matrix_Show(Ai);
        Matrix_Show(Eix);
        Matrix_Show(Eiy);
        for (uint16_t i = 1; i < FitNum - 1; i++)
        {
            Ai->data[i][i - 1] = (dKnotVector[i + 3] * dKnotVector[i + 3]) / (dKnotVector[i + 1] + dKnotVector[i + 2] + dKnotVector[i + 3]);

            Ai->data[i][i] = (dKnotVector[i + 3] * (dKnotVector[i + 1] + dKnotVector[i + 2])) / (dKnotVector[i + 1] + dKnotVector[i + 2] + dKnotVector[i + 3]) +
                    (dKnotVector[i + 2] * (dKnotVector[i + 3] + dKnotVector[i + 4])) / (dKnotVector[i + 2] + dKnotVector[i + 3] + dKnotVector[i + 4]);
            Ai->data[i][i + 1] = (dKnotVector[i + 2] * dKnotVector[i + 2]) / (dKnotVector[i + 2] + dKnotVector[i + 3] + dKnotVector[i + 4]);
            Eix->data[i][0] = (dKnotVector[i + 2] + dKnotVector[i + 3]) * FitVector[i].X;
            Eiy->data[i][0] = (dKnotVector[i + 2] + dKnotVector[i + 3]) * FitVector[i].Y;
        }
        Matrix_Show(Ai);
        Matrix_Show(Eix);
        Matrix_Show(Eiy);
        /* 节点矢量表放大，防止其数据过小导致小于REAL */
        Ai = Matrix_MultConst(Ai, Knot_Max);
        Eix = Matrix_MultConst(Eix, Knot_Max);
        Eiy = Matrix_MultConst(Eiy, Knot_Max);
        Matrix_Show(Ai);
        Matrix_Show(Eix);
        Matrix_Show(Eiy);

        /* 计算控制点 */
        Dix = Matrix_CreateZero(FitNum, 1);
        Diy = Matrix_CreateZero(FitNum, 1);
        AiInverse = Matrix_CreateZero(FitNum, FitNum);
        Matrix_Zero(Dix);
        Matrix_Zero(Diy);
        Matrix_Zero(AiInverse);
        AiInverse = Matrix_Inverse_LU(Ai);
        Matrix_Show(AiInverse);
        Dix = Matrix_Mult(AiInverse, Eix);
        Diy = Matrix_Mult(AiInverse, Eiy);
        Matrix_Show(Dix);
        Matrix_Show(Diy);

        /* 赋值给控制点列表 */
        for (uint16_t i = 0; i < FitNum; i++)
        {
            ctrlPointVector[i + 1].X = Dix->data[i][0];
            ctrlPointVector[i + 1].Y = Diy->data[i][0];
        }

        Matrix_Free(Ai);
        Matrix_Free(Eix);
        Matrix_Free(Eiy);
        Matrix_Free(Dix);
        Matrix_Free(Diy);
        Matrix_Free(AiInverse);

        /* 赋值给控制点列表 */
        ctrlPointVector[0].X = FitVector[0].X;
        ctrlPointVector[0].Y = FitVector[0].Y;
        ctrlPointVector[FitNum + 1].X = FitVector[FitNum - 1].X;
        ctrlPointVector[FitNum + 1].Y = FitVector[FitNum - 1].Y;

        msg = Draw_BsplineWithKnot(level, 0, ratio, ctrlPointVector, ctrlPointNum, KnotVector, KnotNum);
    }
    return msg;
}

/**
 * @brief 绘制B-spline,根据控制点求拟合点
 * @param {REAL} ratio
 * @param {uint16_t} level
 * @param {COORDINATE_t} *ctrlPointVector
 * @param {uint16_t} ctrlPointNum
 * @param {REAL} *KnotVector
 * @param {uint16_t} KnotNum
 * @param {COORDINATE_t} *fitPointVector
 * @param {REAL} speed
 * @return {*}
 * @note
 */
uint16_t MainWindow::Draw_Bspline_GetFitVector(uint16_t level, REAL ratio, REAL speed,
                                               COORDINATE_t *ctrlPointVector, uint16_t ctrlPointNum,
                                               REAL *KnotVector, uint16_t KnotNum,
                                               COORDINATE_t *fitPointVector)
{
    REAL tStep = 0;                  //分辨率
    COORDINATE_t Pxx = {0, 0};       //拟合点
    COORDINATE_t PxxBefore = {0, 0}; //拟合点
    REAL Knot_Max = 0;               //最大节点值
    REAL maxLen = 0;                 //用于机选自适应步长的中间值
    SPLINE_VALUE_t splineMsg = {0, 0};
    QPainter painter(&image); //选入绘图设备中。
    painter.setPen(Qt::yellow);
    REAL len = 0; //长度
    uint16_t fitPointNum = 0; //长度

    SPLINE_MSG_T msg;
    memset(&msg, 0, sizeof(SPLINE_MSG_T));
    /* 判断节点数是否达标 节点数 = 控制点数 + 阶数 + 1 */
    if (KnotNum != ctrlPointNum + level + 1)
    {
        return 0;
    }
    /* 获取最大的节点值 */
    for (uint16_t i = 0; i < KnotNum; ++i)
    {
        if (Knot_Max < KnotVector[i])
        {
            Knot_Max = KnotVector[i];
        }
    }

#if 0
    splineMsg = Draw_Bspline_GetMaxSpeedWithKnot( level,ratio / 100, ctrlPointVector, ctrlPointNum, KnotVector, KnotNum);
    /* 获取最小分辨率的数值的向上取整 */
    maxLen = (int)(splineMsg.value / ratio);
    while (pow(10, tStep) < maxLen)
    {
        tStep++;
    }
    /* 获取最小步长 */
    if (pow(10, tStep) <= maxLen * 2)
    {
        tStep++;
    }
    tStep = 1 / (pow(10, tStep));
    for (REAL t = 0; t < 1; t += tStep)
    {
        Pxx.X = 0;
        Pxx.Y = 0;
        for (uint16_t n = 0; n < ctrlPointNum; n++)
        {
            Pxx.X += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
            Pxx.Y += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
        }
        if (Pxx.X == 0 && Pxx.Y == 0)
        {
            for (uint16_t n = 0; n < ctrlPointNum; n++)
            {
                Pxx.X += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
                Pxx.Y += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
            }
        }

        //        LCD_PutPixel((int16_t)Pxx.X, (int16_t)Pxx.Y);
        if(t != 0)
        {
            len += sqrt((Pxx.X - PxxBefore.X)* (Pxx.X - PxxBefore.X) +(Pxx.Y - PxxBefore.Y)* (Pxx.Y - PxxBefore.Y));
        }
        PxxBefore.X = Pxx.X;
        PxxBefore.Y = Pxx.Y;
        painter.drawPoint(QPoint(Pxx.X, Pxx.Y));
    }
    fitPointNum = (uint16_t)len/speed + 1;//末尾点+1
    tStep = 1 / (REAL)fitPointNum;
    fitPointNum = 0;
#else
    for (uint16_t i = 1; i < ctrlPointNum; i++)
    {
        len += sqrt((ctrlPointVector[i].X - ctrlPointVector[i - 1].X) * (ctrlPointVector[i].X - ctrlPointVector[i - 1].X) +
                (ctrlPointVector[i].Y - ctrlPointVector[i - 1].Y) * (ctrlPointVector[i].Y - ctrlPointVector[i - 1].Y));
    }
#endif
    tStep = speed * ratio / len;

    painter.setPen(Qt::red);

    PxxBefore.X = 0;
    PxxBefore.Y = 0;
    for (REAL t = 0; t < 1; t += tStep)
    {
        Pxx.X = 0;
        Pxx.Y = 0;

        //        for (uint16_t n = 0; n < ctrlPointNum; n++)
        //        {
        //            Pxx.X += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
        //            Pxx.Y += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
        //        }
        //        if (Pxx.X == 0 && Pxx.Y == 0)
        //        {
        //            for (uint16_t n = 0; n < ctrlPointNum; n++)
        //            {
        //                Pxx.X += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
        //                Pxx.Y += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
        //            }
        //        }
        Pxx = Draw_Bspline_DeBoor(level,t,ctrlPointVector,ctrlPointNum,KnotVector,KnotNum);
        if (t != 0)
        {
            len += sqrt((Pxx.X - PxxBefore.X) * (Pxx.X - PxxBefore.X) + (Pxx.Y - PxxBefore.Y) * (Pxx.Y - PxxBefore.Y));
        }
        fitPointVector[fitPointNum++] = Pxx;
        PxxBefore.X = Pxx.X;
        PxxBefore.Y = Pxx.Y;
        painter.drawPoint(QPoint(Pxx.X, Pxx.Y));
        if (t + tStep > 1)
        {
            t = 1;
            Pxx.X = ctrlPointVector[ctrlPointNum - 1].X;
            Pxx.Y = ctrlPointVector[ctrlPointNum - 1].Y;

            // #if !B_SPLINE_DEBOOR
            //             for (uint16_t n = 0; n < ctrlPointNum; n++)
            //             {
            //                 Pxx.X += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
            //                 Pxx.Y += Draw_Bspline_GetBidegN(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
            //             }
            //             if (Pxx.X == 0 && Pxx.Y == 0)
            //             {
            //                 for (uint16_t n = 0; n < ctrlPointNum; n++)
            //                 {
            //                     Pxx.X += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].X;
            //                     Pxx.Y += Draw_Bspline_GetBidegP(level, KnotVector, KnotNum, n, t) * ctrlPointVector[n].Y;
            //                 }
            //             }
            // #else
            //             Pxx = Draw_Bspline_DeBoor(level, t, ctrlPointVector, ctrlPointNum, KnotVector, KnotNum);
            // #endif
            //             if (t != 0)
            //             {
            //                 len += sqrt((Pxx.X - PxxBefore.X) * (Pxx.X - PxxBefore.X) + (Pxx.Y - PxxBefore.Y) * (Pxx.Y - PxxBefore.Y));
            //             }
            //             PxxBefore.X = Pxx.X;
            //             PxxBefore.Y = Pxx.Y;
            fitPointVector[fitPointNum++] = Pxx;
            painter.drawPoint(QPoint(Pxx.X, Pxx.Y));
            break;
        }
    }

    /* 显示最小矩形 */
    if (B_SPLINE_SHOW_CTRL)
    {
    }
    painter.end();
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
    return fitPointNum;
}

void MainWindow::on_pushButton_clear_clicked()
{
    ctrlPointCnt = 0;
    memset(ctrlPointForBspline, 0, sizeof(ctrlPointForBspline));
    ui->comboBox_lineType->setEnabled(true);
    image.fill(Qt::white);
    my_scene->addPixmap(image);
    ui->canvas->setScene(my_scene);
    ui->canvas->show();
    ui->textBrowser_X->clear();
    ui->textBrowser_Y->clear();
    ui->textBrowser_msg->clear();
    ui->textBrowser_X->show();
    ui->textBrowser_Y->show();
    ui->textBrowser_msg->show();
}
