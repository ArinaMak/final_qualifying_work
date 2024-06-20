#ifndef TASK_H
#define TASK_H
#include <QtMath>
#include <QDebug>

// #include "mainwindow.h"

class Task
{
public:
    double l = 3.0;             // Длина стержня
    double F_2=20.0;
    double F_1=10.0;
    double M = 0.0;         // Сосредоточенная масса
    double omega = 5299.9;      // Частота колебаний

    double rho = 7810.0;        // Плотность материала стержня
    double E = 200e9;           // Модуль Юнга

    double a;                   // a=qSqrt(E/ro)
    double beta;                // beta=omega/a

    double delta = 0.01;
    double eps=0.001;

    int N=100;

    int k_min=-1;
    int k_max;

    int k_curr=k_min;
    QVector<double> x_1_curr_interval;
    QVector<double> x_2_curr_interval;

    QVector<double> x_1_arr;
    QVector<double> x_2_arr;

    double x_opt1=0.0;
    double x_opt2=0.0;

    int task_number=0;
    int is_task_success=0;

    Task(double il, double iF_1, double iF_2, double iM, double iomega, double irho, double iE, double ieps, int iN, int itask_number);
    double roundingNumber(double x, int n);

    QVector<double> der_OnePart(double x);

    bool isTheOptimalPoint_OnePart(double x);

    void new_M_and_l_FromSpecifiedF_OnePart();
    void new_l_and_F2_FromSpecifiedF1_and_M_OnePart();
    void new_l_and_F1_FromSpecifiedF2_and_M_OnePart();
    void new_M_and_F1_FromSpecifiedF2_and_l_OnePart();
    void new_M_and_F2_FromSpecifiedF1_and_l_OnePart();

    void fillingInVectors_OnePart(QVector<double> &x_i, QVector<double> &F_i);

    // a rod consisting of two parts - a variable and a constant
    double func_l_TwoParts1(double x_1, double x_l);

    double searchFor_x_TwoParts1();
    bool isTheOptimalPoint_TwoParts1(double x);
    QVector<double> decision_TwoParts1();
    void fillingInVectors_TwoParts1(QVector<double> &x_i, QVector<double> &F_i, double x_2);

    // double func_l_TwoParts1(double x_1, double x_l);
    void findOptimalLength_TwoParts1(double x_2);


    // a rod consisting of two parts - a constant and a variable
    double f_TwoParts2(double x, double alpha);
    double searchFor_x_usingTheHalfDivisionMethodOnTheCurrentInterval_TwoParts2(double left_border, double right_border);
    bool isTheOptimalPoint_TwoParts2(double x);

    QVector<double> decision_TwoParts2();
    void fillingInVectors_TwoParts2(QVector<double> &x_i, QVector<double> &F_i, double x_1);

    void findOptimalLength_TwoParts2();



    void fillingInVectors(QVector<double> &x_i, QVector<double> &F_i);


    QVector<double> f_k(double x_1, double x_2);
    QVector<double> matrixJacobi(double x_1, double x_2);
    QVector<double> inverseMatrixJacobi(QVector<double> J);
    QVector<double> matrixVectorMultiplication(QVector<double> matrixA, int colsA, int rowsA, QVector<double> matrixB, int colsB, int rowsB);

    double findingAlpha(double x_1);
    double findingAlpha2(double x_2);
    QVector<double> decisionThreeParts();
    // QVector<double> decisionThreeParts2(Ui::MainWindow *ui);
    QVector<double> methodNewton(QVector<double> x_k);

    void fillingInVectors3(QVector<double> &x_i, QVector<double> &F_i, QVector<double> x_opt);



    int searching_n_ForTheInterval_x_2(int k);
    QVector<double> searching_Interval_x_1(int k);
    QVector<double> searching_Interval_x_2(int n);


    QVector<double> findingInitialApproximationFromTheCorner(double corner, bool is_left_root);

    double Norm(QVector<double> x_i, QVector<double> x_i_prev);


    bool checkingTheResultsWithZeroMassAndTheTaskFromTheBook(double x_1, double x_2, double alpha);
};

#endif // TASK_H
