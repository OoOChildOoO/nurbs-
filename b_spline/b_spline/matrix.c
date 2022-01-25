#include "matrix.h"

#if 0
// 初始化矩阵(将所有元素初始化为０)
void Matrix_Zero(Matrix mat)
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

void Matrix_Set(Matrix mat, int row, int col, REAL value)
{
    if (row < mat->row && col < mat->column)
        mat->data[row][col] = value;
}

//　创建矩阵
Matrix Matrix_CreateZero(int row, int col)
{
    Matrix mat;
    mat = (Matrix)malloc(sizeof(struct MNode)); //　分配结构体指针
    if (row <= 0 || col <= 0)
    {
        //        printf("ERROR, in creat_Matrix the row or col <= 0\n");
        exit(1);
    }
    if (row > 0 && col > 0)
    {
        mat->row = row;
        mat->column = col;
        mat->data = (REAL **)malloc(row * sizeof(REAL *)); // 分配头指针
        if (mat->data == NULL)
        {
            //            printf("ERROR, in creat_Matrix the mat->data == NULL\n");
            exit(1);
        }
        int i;
        for (i = 0; i < row; i++)
        {
            *(mat->data + i) = (REAL *)malloc(col * sizeof(REAL)); //　分配每行的指针
            if (mat->data[i] == NULL)
            {
                //                printf("ERROR, in create_Matrix the mat->data[i] == NULL\n");
                exit(1);
            }
        }
        Matrix_Zero(mat);
    }
    return mat;
}

// 释放申请的矩阵空间
void Matrix_Free(Matrix mat)
{
    for (int i = 0; i < mat->row; i++)
        free(mat->data[i]); // 释放行指针
    free(mat->data);        // 释放头指针
    free(mat);              // 释放结构体指针
}

//　创建单位矩阵
Matrix Matrix_CreateEye(int n)
{
    Matrix E;
    int i, j;
    if (n <= 0)
    {
        //        printf("ERROR in Matrix_CreateEye\n");
        exit(1);
    }
    E = Matrix_CreateZero(n, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
                E->data[i][j] = 1;
            else
                E->data[i][j] = 0;
        }
    }
    return E;
}

//　打印矩阵
void Matrix_Show(Matrix mat)
{
    int i, j;
    //    printf("%s\n", s);
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
        {
            //            printf("%.6f\t", mat->data[i][j]);
            qDebug() << mat->data[i][j];
        }
        //        printf("\n");
        qDebug() << "\n";
    }
    //    printf("\n");
    qDebug() << "\n";
}

// 给矩阵每个元素赋值
void Matrix_SetData(Matrix mat, REAL data[])
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

//　取出矩阵某行某列的元素
REAL Matrix_At(Matrix mat, int row, int col)
{
    REAL rst;
    rst = mat->data[row - 1][col - 1];
    return rst;
}

// 矩阵第n行与第m行互换
void Matrix_Swap(Matrix mat, int n, int m)
{
    REAL temp;
    for (int i = 0; i < mat->column; i++)
    {
        temp = mat->data[n - 1][i];
        mat->data[n - 1][i] = mat->data[m - 1][i];
        mat->data[m - 1][i] = temp;
    }
}

//　对一个矩阵进行复制
Matrix Matrix_Copy(Matrix mat)
{
    Matrix copy_mat = Matrix_CreateZero(mat->row, mat->column);
    for (int i = 0; i < mat->row; i++)
    {
        for (int j = 0; j < mat->column; j++)
            copy_mat->data[i][j] = mat->data[i][j];
    }
    return copy_mat;
}

//　伴随矩阵
Matrix Matrix_Adjoint(Matrix mat)
{
    //    Matrix adj;
    return mat;
}

//　矩阵加法
Matrix Matrix_Add(Matrix mat_1, Matrix mat_2)
{
    Matrix rst_mat;
    if (mat_1->column != mat_2->column)
    {
        //        printf("ERROR in AddorSub, column !=\n");
        exit(1);
    }
    if (mat_1->row != mat_2->row)
    {
        printf("ERROR in AddorSub, row !=\n");
        exit(1);
    }
    int i, j;
    rst_mat = Matrix_CreateZero(mat_1->row, mat_1->column);
    for (i = 0; i < mat_1->row; i++)
    {
        for (j = 0; j < mat_1->column; j++)
            rst_mat->data[i][j] = mat_1->data[i][j] + mat_2->data[i][j];
    }
    return rst_mat;
}

//  矩阵的减法
Matrix Matrix_Sub(Matrix mat_1, Matrix mat_2)
{
    Matrix rst_mat;
    if (mat_1->column != mat_2->column)
    {
        //        printf("ERROR in AddorSub, column !=\n");
        exit(1);
    }
    if (mat_1->row != mat_2->row)
    {
        printf("ERROR in AddorSub, row !=\n");
        exit(1);
    }
    int i, j;
    rst_mat = Matrix_CreateZero(mat_1->row, mat_1->column);
    for (i = 0; i < mat_1->row; i++)
    {
        for (j = 0; j < mat_1->column; j++)
            rst_mat->data[i][j] = mat_1->data[i][j] - mat_2->data[i][j];
    }
    return rst_mat;
}

// 矩阵转置
Matrix Matrix_Transpose(Matrix mat)
{
    Matrix mat_;
    int i, j;
    mat_ = Matrix_CreateZero(mat->row, mat->column);
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
            mat_->data[i][j] = mat->data[i][j];
    }
    return mat_;
}

// 矩阵乘法
Matrix Matrix_Mult(Matrix mat_1, Matrix mat_2)
{
    Matrix rst_mat;
    int i, j, m;
    if (mat_1->column != mat_2->row)
    {
        //        printf("ERROR in Matrix_Mult, column != row\n");
        exit(1);
    }
    else
    {
        rst_mat = Matrix_CreateZero(mat_1->row, mat_2->column);
        for (i = 0; i < mat_1->row; i++)
        {
            for (j = 0; j < mat_2->column; j++)
            {
                for (m = 0; m < mat_1->column; m++)
                    rst_mat->data[i][j] += mat_1->data[i][m] * mat_2->data[m][j];
            }
        }
    }
    return rst_mat;
}

// 矩阵乘法常数
Matrix Matrix_MultConst(Matrix mat, float ratio)
{
    int i = 0, j = 0;
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
        {
            mat->data[i][j] *= ratio;
        }
    }
    return mat;
}

// 矩阵求逆，利用矩阵LU分解求逆（直接三角分解）
Matrix Matrix_Inverse_LU(Matrix mat)
{
    Matrix inv;
    if (mat->column != mat->row)
    {
        //        printf("ERROR in inverse, the row != the column\n");
        exit(1);
    }
    if (Matrix_Det(mat) == 0)
    {
        //        printf("The Matrix is not invertible\n");
        exit(1);
    }
    // Ｌ矩阵和Ｕ矩阵
    Matrix L, U;
    int n = mat->column;
    inv = Matrix_CreateZero(n, n);
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
    Matrix L_, U_;
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
    inv = Matrix_Mult(U_, L_);

    Matrix_Free(L_);
    Matrix_Free(U_);
    Matrix_Free(L);
    Matrix_Free(U);
    return inv;
}

// 矩阵求逆，利用初等行变换求逆
Matrix Matrix_Inverse_EleTrans(Matrix mat)
{
    if (mat->column != mat->row)
    {
        //        printf("ERROR in inverse, the row != the column\n");
        exit(1);
    }
    if (Matrix_Det(mat) == 0)
    {
        //        printf("The Matrix is not invertible\n");
        exit(1);
    }
    int n = mat->row;
    // 为防止改变原矩阵，此处进行复制处理
    Matrix mat_ = Matrix_Copy(mat);
    //　创建单位矩阵
    Matrix inv_eye = Matrix_CreateEye(n);
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
                temp_e = inv_eye->data[k][j];
                inv_eye->data[k][j] = inv_eye->data[cnt][j];
                inv_eye->data[cnt][j] = temp_e;
            }
        }
        // 消元
        for (i = k + 1; i < n; i++)
        {
            e = mat_->data[i][k] / mat_->data[k][k];
            for (j = 0; j < n; j++)
            {
                mat_->data[i][j] = mat_->data[i][j] - e * mat_->data[k][j];
                inv_eye->data[i][j] = inv_eye->data[i][j] - e * inv_eye->data[k][j];
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
            inv_eye->data[i][j] = inv_eye->data[i][j] * e;
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
                inv_eye->data[j][t] = inv_eye->data[j][t] - e * inv_eye->data[i][t];
            }
        }
    }
    Matrix_Free(mat_);
    return inv_eye;
}

//　顺序高斯消去法
Matrix Matrix_Elimination_GaussOrdinal(Matrix A, Matrix b)
{
    int n;
    n = b->row;
    Matrix x;
    x = Matrix_CreateZero(b->row, b->column);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            if (A->data[i][j] == 0)
            {
                //                printf("can't use the Matrix_Elimination_GaussOrdinal\n");
                exit(1);
            }
    }
    // 消元
    REAL L[n];
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            L[i] = A->data[i][k] / A->data[k][k];
            for (int j = k + 1; j < n; j++)
                A->data[i][j] = A->data[i][j] - L[i] * A->data[k][j];
            b->data[i][0] = b->data[i][0] - L[i] * b->data[k][0];
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (A->data[i][i] == 0)
        {
            //            printf("can't use the Matrix_Elimination_GaussOrdinal\n");
            exit(1);
        }
    }
    // 回代
    x->data[n - 1][0] = b->data[n - 1][0] / A->data[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        REAL sum_a = 0;
        for (int j = i + 1; j < n; j++)
            sum_a += A->data[i][j] * x->data[j][0];
        x->data[i][0] = (b->data[i][0] - sum_a) / A->data[i][i];
    }
    return x;
}

//　列主元高斯消去法
Matrix Matrix_Elimination_GaussColumn(Matrix A, Matrix b)
{
    int n;
    n = b->row;
    Matrix x;
    x = Matrix_CreateZero(b->row, b->column);
    REAL L[n];
    for (int k = 0; k < n - 1; k++)
    {
        int cnt = k;
        REAL a_max = fabs(A->data[k][k]);
        for (int i = k; i < n; i++)
        {
            if (fabs(A->data[i][k]) > a_max)
            {
                a_max = fabs(A->data[i][k]);
                cnt = i; // 确定下标i
            }
        }
        if (A->data[cnt][k] == 0)
        {
            //            printf("Matrix_Elimination_GaussColumn: no unique solution\n");
            exit(1);
        }
        // 换行
        if (cnt != k)
        {
            REAL t = 0, s = 0;
            for (int j = k; j < n; j++)
            {
                t = A->data[k][j];
                A->data[k][j] = A->data[cnt][j];
                A->data[cnt][j] = t;
                s = b->data[cnt][0];
                b->data[cnt][0] = b->data[k][0];
                b->data[k][0] = s;
            }
        }
        // step 5
        for (int i = k + 1; i < n; i++)
        {
            L[i] = A->data[i][k] / A->data[k][k];
            for (int j = k + 1; j < n; j++)
                A->data[i][j] = A->data[i][j] - L[i] * A->data[k][j];
            b->data[i][0] = b->data[i][0] - L[i] * b->data[k][0];
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (A->data[i][i] == 0.0)
        {
            //            printf("Matrix_Elimination_GaussColumn: no unique solution\n");
            exit(1);
        }
    }
    // 回代
    x->data[n - 1][0] = b->data[n - 1][0] / A->data[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        REAL sum_a = 0;
        for (int j = i + 1; j < n; j++)
            sum_a += A->data[i][j] * x->data[j][0];
        x->data[i][0] = (b->data[i][0] - sum_a) / A->data[i][i];
    }
    return x;
}

// 通过行变换（列主元消去法）将矩阵变换成上三角矩阵.注意行变换会改变行列式的符号
REAL Matrix_Det(Matrix mat)
{
    REAL det = 1.0;
    if (mat->row != mat->column)
    {
        //        printf("ERROR in Matrix_Det, the row != the column\n");
        exit(1);
    }
    if (mat->row == 1 && mat->column == 1)
        return mat->data[0][0];
    // 为防止改变原矩阵，此处进行复制处理
    Matrix mat_ = Matrix_Copy(mat);
    int n = mat_->row, s = 0, cnt;
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

#else
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

// /**
//  * @brief 顺序高斯消去法
//  * @param {_IN Matrix} mat_1
//  * @param {_IN Matrix} mat_2
//  * @return {*}
//  * @note
//  */
// Matrix Matrix_Elimination_GaussOrdinal(_IN Matrix mat_1, _IN Matrix mat_2)
// {
//     MATRIX_RES_T res = MATRIX_NO_ERROR;
//     Matrix matRes;
//     int n;
//     REAL L[mat_2->row];
//     n = mat_2->row;
//     matRes = Matrix_CreateZero(mat_2->row, mat_2->column);
//     for (int i = 0; i < n; i++)
//     {
//         for (int j = 0; j < n; j++)
//             if (mat_1->data[i][j] == 0)
//             {
//                 res = MATRIX_NOT_GAUSSORDINAL;
//                 return NULL;
//             }
//     }
//     // 消元
//     for (int k = 0; k < n - 1; k++)
//     {
//         for (int i = k + 1; i < n; i++)
//         {
//             L[i] = mat_1->data[i][k] / mat_1->data[k][k];
//             for (int j = k + 1; j < n; j++)
//                 mat_1->data[i][j] = mat_1->data[i][j] - L[i] * mat_1->data[k][j];
//             mat_2->data[i][0] = mat_2->data[i][0] - L[i] * mat_2->data[k][0];
//         }
//     }
//     for (int i = 0; i < n; i++)
//     {
//         if (mat_1->data[i][i] == 0)
//         {
//             res = MATRIX_NOT_GAUSSORDINAL;
//             return NULL;
//         }
//     }
//     // 回代
//     matRes->data[n - 1][0] = mat_2->data[n - 1][0] / mat_1->data[n - 1][n - 1];
//     for (int i = n - 2; i >= 0; i--)
//     {
//         REAL sum_a = 0;
//         for (int j = i + 1; j < n; j++)
//             sum_a += mat_1->data[i][j] * matRes->data[j][0];
//         matRes->data[i][0] = (mat_2->data[i][0] - sum_a) / mat_1->data[i][i];
//     }
//     return matRes;
// }

// /**
//  * @brief 列主元高斯消去法
//  * @param {_IN Matrix} mat_1
//  * @param {_IN Matrix} mat_2
//  * @return {*}
//  * @note 
//  */
// Matrix Matrix_Elimination_GaussColumn(_IN Matrix mat_1, _IN Matrix mat_2)
// {
//     MATRIX_RES_T res = MATRIX_NO_ERROR;
//     Matrix matRes;
//     int n;
//     REAL L[mat_2->row];
//     n = mat_2->row;
//     matRes = Matrix_CreateZero(mat_2->row, mat_2->column);
//     for (int k = 0; k < n - 1; k++)
//     {
//         int cnt = k;
//         REAL a_max = fabs(mat_1->data[k][k]);
//         for (int i = k; i < n; i++)
//         {
//             if (fabs(mat_1->data[i][k]) > a_max)
//             {
//                 a_max = fabs(mat_1->data[i][k]);
//                 cnt = i; // 确定下标i
//             }
//         }
//         if (mat_1->data[cnt][k] == 0)
//         {
//             res = MATRIX_NOT_GAUSSORDINAL;
//             return NULL;
//         }
//         // 换行
//         if (cnt != k)
//         {
//             REAL t = 0, s = 0;
//             for (int j = k; j < n; j++)
//             {
//                 t = mat_1->data[k][j];
//                 mat_1->data[k][j] = mat_1->data[cnt][j];
//                 mat_1->data[cnt][j] = t;
//                 s = mat_2->data[cnt][0];
//                 mat_2->data[cnt][0] = mat_2->data[k][0];
//                 mat_2->data[k][0] = s;
//             }
//         }
//         // step 5
//         for (int i = k + 1; i < n; i++)
//         {
//             L[i] = mat_1->data[i][k] / mat_1->data[k][k];
//             for (int j = k + 1; j < n; j++)
//                 mat_1->data[i][j] = mat_1->data[i][j] - L[i] * mat_1->data[k][j];
//             mat_2->data[i][0] = mat_2->data[i][0] - L[i] * mat_2->data[k][0];
//         }
//     }
//     for (int i = 0; i < n; i++)
//     {
//         if (mat_1->data[i][i] == 0.0)
//         {
//             res = MATRIX_NOT_GAUSSORDINAL;
//             return NULL;
//         }
//     }
//     // 回代
//     matRes->data[n - 1][0] = mat_2->data[n - 1][0] / mat_1->data[n - 1][n - 1];
//     for (int i = n - 2; i >= 0; i--)
//     {
//         REAL sum_a = 0;
//         for (int j = i + 1; j < n; j++)
//             sum_a += mat_1->data[i][j] * matRes->data[j][0];
//         matRes->data[i][0] = (mat_2->data[i][0] - sum_a) / mat_1->data[i][i];
//     }
//     return matRes;
// }

#endif
