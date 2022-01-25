/*
 * @Author your name
 * @Date 2021-12-17 11:00:33
 * @LastEditTime 2021-12-20 16:56:40
 * @LastEditors Please set LastEditors
 * @Description 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath \undefinede:\WJ\JK\Program\T3020M\panel\learn\b_spline_learn\b_spline\b_spline\matrix.h
 */
#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define _IN
#define _OUT

typedef float REAL;
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


#endif // MATRIX_H
