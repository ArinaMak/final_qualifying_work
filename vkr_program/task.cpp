#include "task.h"
#include <QtMath>

#define M_PI           3.141592653589793238462643383279502884197  /* pi */

Task::Task(double il, double iF_1, double iF_2, double iM, double iomega, double irho, double iE, double ieps, int iN, int itask_number)
{
    l = il;             // Длина стержня
    F_1=iF_1;
    F_2=iF_2;
    M=iM;
    omega=iomega;
    rho=irho;
    E=iE;
    eps=ieps;
    eps=0.0001;
    N=iN;
    task_number=itask_number;

    a=qSqrt(E/rho);
    beta=omega/a;
    k_max=floor(l*beta/M_PI-3.0/4.0);
    k_curr=k_min;

    int n=searching_n_ForTheInterval_x_2(k_curr);
    x_1_curr_interval = searching_Interval_x_1(k_curr);
    x_2_curr_interval = searching_Interval_x_2(n);


    qDebug()<<"k_min = "<<k_min<<"k_max = "<<k_max;
    qDebug()<<searching_Interval_x_1(k_min);
    qDebug()<<searching_Interval_x_1(k_max);
    qDebug()<<"l = "<<l<<"\nF_1 = "<<F_1<<"\nF_2 = "<<F_2<<"\nM = "<<M<<"\nomega = "<<omega<<"\nrho = "<<rho<<"\nE = "<<E<<"\na = "<<a<<"beta = "<<beta;
}

double Task::roundingNumber(double x, int n)
{
    QString tmp_x_str=QString("%1").arg(x, 0, 'f', n);
    return tmp_x_str.toDouble();
}


double Task::func_l_TwoParts1(double x_1, double x_l)
{
    double tmp=tanh(beta*x_1)-1/(tan(beta*(x_l-x_1)+atan(M*beta*a*a/(E*F_1))));
    return tmp;
}


QVector<double> Task::der_OnePart(double x)
{
    QVector<double> der_tmp;

    der_tmp.push_back(-M*beta*beta*a*a/(cosh(beta*x)*cosh(beta*x)));
    der_tmp.push_back(F_1*2.0*beta*cosh(beta*x)*sinh(beta*x));
    return der_tmp;
}



bool Task::isTheOptimalPoint_OnePart(double x)
{
    // qDebug()<<QString("M = %1, \nl = %2, \nF_1 = %3, \nF_2 = %4").arg(l, 0, 'f', 16).arg(M, 0, 'f', 16).arg(F_1, 0, 'f', 16).arg(F_2, 0, 'f', 16);

    double tmp = M*beta*a*a*tanh(beta*x)-E;
    double tmp2 = F_1*(cosh(beta*x)*cosh(beta*x))-F_2;

    //qDebug()<<"true123\neps = "+QString("%1, \nfabs(tmp) = %2,\nfabs(tmp2) = %3").arg(eps, 0, 'f', 16).arg(fabs(tmp), 0, 'f', 16).arg(fabs(tmp2), 0, 'f', 16);
    QVector<double> der_tmp = der_OnePart(x);
    qDebug()<<"der_tmp1 = "+QString("%1, \nder_tmp2 = %23").arg(der_tmp.value(0), 0, 'f', 16).arg(der_tmp.value(1), 0, 'f', 16);
    qDebug()<<"true123\neps = "+QString("%1, \nfabs(tmp) = %2,\nfabs(tmp2) = %3").arg(eps, 0, 'f', 16).arg(fabs(tmp), 0, 'f', 16).arg(fabs(tmp2), 0, 'f', 16);

    double eps1=1e-16;
    qDebug()<<"der_tmp123 = "+QString("%1, \nder_tmp2345 = %23").arg(der_tmp.value(0)*eps1, 0, 'f', 16).arg(der_tmp.value(1)*eps1, 0, 'f', 16);
    if(fabs(tmp)<eps && fabs(tmp2)<eps) return true;
    // qDebug()<<"isTheOptimalPoint_OnePart(double x) = false";
    return false;
}

void Task::new_M_and_l_FromSpecifiedF_OnePart()
{
    l=acosh(sqrt(F_2/F_1))/beta;
    M=E/(beta*a*a*tanh(l*beta));
}

void Task::new_l_and_F2_FromSpecifiedF1_and_M_OnePart()
{
    if(M==0 || E/(M*beta*a*a)>=1)
    {
        l=0.0;
        return;
    }
    l=atanh(E/(M*beta*a*a))/beta;
    F_2=F_1*(cosh(beta*l)*cosh(beta*l));
}

void Task::new_l_and_F1_FromSpecifiedF2_and_M_OnePart()
{
    if(M==0 || E/(M*beta*a*a)>=1)
    {
        qDebug()<<"!!!\nE/(M*beta*a*a) = "<<E/(M*beta*a*a);
        l=0.0;
        return;
    }
    l=atanh(E/(M*beta*a*a))/beta;
    F_1=F_2/(cosh(beta*l)*cosh(beta*l));
}

void Task::new_M_and_F1_FromSpecifiedF2_and_l_OnePart()
{
    M=E/(beta*a*a*tanh(beta*l));
    F_1=F_2/(cosh(beta*l)*cosh(beta*l));
}

void Task::new_M_and_F2_FromSpecifiedF1_and_l_OnePart()
{
    M=E/(beta*a*a*tanh(beta*l));
    F_2=F_1*(cosh(beta*l)*cosh(beta*l));
}

void Task::fillingInVectors_OnePart(QVector<double> &x_i, QVector<double> &F_i)
{
    double h=(l-0)/(N-1);

    for(int i=0; i<N;i++)
    {
        x_i.push_back(0+i*h);

        double F_tmp=F_2/(cosh(beta*x_i.value(i))*cosh(beta*x_i.value(i)));
        F_i.push_back(F_tmp);

    }
}

double Task::searchFor_x_TwoParts1()
{
    double tmp_x =roundingNumber(acosh(sqrt(F_2/F_1))/beta, 3);
    return tmp_x;
    // return acosh(sqrt(F_2/F_1))/beta;
}

bool Task::isTheOptimalPoint_TwoParts1(double x)
{
    double tmp0=atan(M*beta*a*a/(E*F_1));
    double tmp=tanh(beta*x)-1.0/(tan(beta*(l-x)+tmp0));

    if(fabs(tmp-0.0)<delta && 0.0<x) return true;
    return false;
}

QVector<double> Task::decision_TwoParts1()
{
    QVector<double> x_2;

    double x_2_tmp=searchFor_x_TwoParts1();
    QString tmp_x_2_str=QString("%1").arg(x_2_tmp, 0, 'f', 3);
    x_2_tmp=tmp_x_2_str.toDouble();

    if(isTheOptimalPoint_TwoParts1(x_2_tmp)) x_2.push_back(x_2_tmp);

    return x_2;
}

void Task::fillingInVectors_TwoParts1(QVector<double> &x_i, QVector<double> &F_i, double x_2)
{
    double h=(x_2-0)/(N-2);

    for(int i=0; i<N-1;i++)
    {
        x_i.push_back(0+i*h);
        double F_tmp=F_2/(cosh(beta*x_i.value(i))*cosh(beta*x_i.value(i)));
        F_i.push_back(F_tmp);
    }
    x_i.push_back(l);
    F_i.push_back(F_1);
}

void Task::findOptimalLength_TwoParts1(double x_2)
{
    // double
    int k = ceil( atan(M*beta*a*a/(E*F_1))/M_PI );

    double left = 0.0;
    double right = 0.0;

    if(k!=atan(M*beta*a*a/(E*F_1))/M_PI)
    {
        left = x_2;
        right = ( M_PI*k-atan(M*beta*a*a/(E*F_1))+beta*x_2 )/beta;
    }
    else
    {
        left=( M_PI*k-atan(M*beta*a*a/(E*F_1))+beta*x_2 )/beta;
        right = ( M_PI*(k+1)-atan(M*beta*a*a/(E*F_1))+beta*x_2 )/beta;
    }

    // qDebug()<<QString("atan(M*beta*a*a/(E*F_1))/M_PI=%1\nceil( atan(M*beta*a*a/(E*F_1))/M_PI ) = %2").arg(atan(M*beta*a*a/(E*F_1))/M_PI, 0, 'f', 16).arg(ceil( atan(M*beta*a*a/(E*F_1))/M_PI ), 0, 'f', 16);
    // qDebug()<<QString("x_2=%1\nright = %2").arg(left, 0, 'f', 16).arg(right, 0, 'f', 16);

    double left_border=left;
    double right_border=right;

    left=left+eps;
    right=right-eps;

    double x_0=0.0;

    if((func_l_TwoParts1(x_2, left)*func_l_TwoParts1(x_2, right)) > 0)
    {
        qDebug()<<"На заданном интервале нет корней";
        return;
    }

    x_0=(left+right)/2.0;
    while(fabs(left-right)>=eps)
    {
        if((func_l_TwoParts1(x_2, left)*func_l_TwoParts1(x_2, x_0))>0) left=x_0;
        else right=x_0;

        x_0=(left+right)/2.0;
    }
    if(left_border<x_0 && x_0<right_border) l= roundingNumber(x_0, 3);
    else l=0.0;
    if(!isTheOptimalPoint_TwoParts1(x_2)) l=0.0;
}

double Task::f_TwoParts2(double x, double alpha)
{
    return tan(beta*x)-tanh(beta*x+alpha);
}

QVector<double> Task::decision_TwoParts2()
{
    k_max=floor(l*beta/M_PI-1.0/2.0);
    k_curr=k_min;

    double x_1_tmp=0.0;
    QVector<double> x_1;

    if(k_max<k_min)
    {
        x_1_tmp = searchFor_x_usingTheHalfDivisionMethodOnTheCurrentInterval_TwoParts2(0, l);
        x_1_tmp=roundingNumber(x_1_tmp, 3);
        if(isTheOptimalPoint_TwoParts2(x_1_tmp)) x_1.push_back(x_1_tmp);
    }
    else
    {
        double left=0.0;
        double right=(M_PI/2.0+M_PI*k_min)/beta;
        x_1_tmp = searchFor_x_usingTheHalfDivisionMethodOnTheCurrentInterval_TwoParts2(left, right);
        x_1_tmp=roundingNumber(x_1_tmp, 3);

        if(isTheOptimalPoint_TwoParts2(x_1_tmp)) x_1.push_back(x_1_tmp);
        while(k_curr<k_max)
        {
            left=(M_PI/2.0+M_PI*k_curr)/beta;
            right=(M_PI/2.0+M_PI*(k_curr+1))/beta;
            x_1_tmp = searchFor_x_usingTheHalfDivisionMethodOnTheCurrentInterval_TwoParts2( left, right);

            if(isTheOptimalPoint_TwoParts2(x_1_tmp)) x_1.push_back(x_1_tmp);
            k_curr++;
        }
        left=(M_PI/2.0+M_PI*k_max)/beta;
        right=l;
        x_1_tmp = searchFor_x_usingTheHalfDivisionMethodOnTheCurrentInterval_TwoParts2(left, right);
        x_1_tmp=roundingNumber(x_1_tmp, 3);

        if(isTheOptimalPoint_TwoParts2(x_1_tmp)) x_1.push_back(x_1_tmp);
    }

    return x_1;
}

bool Task::isTheOptimalPoint_TwoParts2(double x)
{
    double alpha=atanh(E*F_1/(M*beta*a*a))-beta*l;
    double tmp = F_2*(cosh(beta*x+alpha)*cosh(beta*x+alpha))-F_1*(cosh(beta*l+alpha)*cosh(beta*l+alpha));
    // qDebug()<<"\nfabs(tmp-0.0) = "<<fabs(tmp-0.0);
    qDebug()<<"SUCCESS??? tmp = "<<tmp;
    if(fabs(tmp-0.0)<delta && 0.0<x) return true;
    return false;
}

void Task::fillingInVectors_TwoParts2(QVector<double> &x_i, QVector<double> &F_i, double x_1)
{
    double h=(l-x_1)/(N-2);
    double alpha=atanh(E*F_1/(M*beta*a*a))-beta*l;
    double D=F_1*cosh(beta*l+alpha)*cosh(beta*l+alpha);

    x_i.push_back(0);
    F_i.push_back(F_2);
    for(int i=0; i<N-1;i++)
    {
        x_i.push_back(x_1+i*h);
        double F_tmp=D/(cosh(beta*x_i.value(i+1)+alpha)*cosh(beta*x_i.value(i+1)+alpha));
        F_i.push_back(F_tmp);
    }
}

void Task::findOptimalLength_TwoParts2()
{
    qDebug()<<"0) E*F_1/(M*beta*a*a) = "<<E*F_1/(M*beta*a*a);
    double tmp = -F_2*(1-(tanh(atanh(E*F_1/(M*beta*a*a))))*(tanh(atanh(E*F_1/(M*beta*a*a)))))/F_1+1.0;
    qDebug()<<"1) tmp = "<<tmp;
    if(tmp<0) return;
    tmp=sqrt(tmp);
    qDebug()<<"2) sqrt(tmp) = "<<tmp;
    //   //tmp=atan(tmp)/beta;
    tmp=atan(tmp)/beta;
    // tmp=(atan(tmp)+M_PI)/beta;
    double tmp2=tmp+M_PI/beta;
    double l_tmp2=0.0;

    tmp=roundingNumber(tmp, 3);
    tmp2=roundingNumber(tmp2, 3);
    double alpha=atanh(tan(beta*tmp))-beta*tmp;
    double l_tmp = ( atanh(E*F_1/(M*beta*a*a))-alpha )/beta;
    alpha=atanh(tan(beta*tmp2))-beta*tmp2;
    l_tmp2 = ( atanh(E*F_1/(M*beta*a*a))-alpha )/beta;
    while(l_tmp2<l)
    {
        tmp=tmp2;
        l_tmp=l_tmp2;
        qDebug()<<"!!! tmp2 = "<<tmp2;
        tmp2+=M_PI/beta;
        tmp2=roundingNumber(tmp2, 3);
        qDebug()<<"!!! tmp2 = "<<tmp2;
        alpha=atanh(tan(beta*tmp2))-beta*tmp2;
        l_tmp2 = ( atanh(E*F_1/(M*beta*a*a))-alpha )/beta;
    }
    qDebug()<<"l_tmp = "<<l_tmp<<"\n\nl_tmp2 = "<<l_tmp2;
    if(fabs(l_tmp-l)>fabs(l_tmp2-l)) x_opt1=tmp2;
    else x_opt1=tmp;
    qDebug()<<"!!! x_opt1 = "<<x_opt1;
    // qDebug()<<"3) atan(tmp)/beta = "<<tmp;
    // x_opt1=roundingNumber(tmp, 3);

    alpha=atanh(tan(beta*x_opt1))-beta*x_opt1;
    l_tmp = ( atanh(E*F_1/(M*beta*a*a))-alpha )/beta;

    l=roundingNumber(l_tmp, 3);

    if(!isTheOptimalPoint_TwoParts2(x_opt1)) l=0.0;
}

double Task::searchFor_x_usingTheHalfDivisionMethodOnTheCurrentInterval_TwoParts2(double left_border, double right_border)
{
    double left=left_border+eps;
    double right=right_border-eps;
    double alpha=atanh(E*F_1/(M*beta*a*a))-beta*l;
    double x_0=0.0;

    if((f_TwoParts2(left, alpha)*f_TwoParts2(right, alpha)) > 0) return 0.0;

    x_0=(left+right)/2.0;
    while(fabs(left-right)>=eps)
    {
        if((f_TwoParts2(left, alpha)*f_TwoParts2(x_0, alpha))>0) left=x_0;
        else right=x_0;

        x_0=(left+right)/2.0;
    }
    if(left_border<x_0 && x_0<right_border) return x_0;
    return 0.0;
}



void Task::fillingInVectors(QVector<double> &x_i, QVector<double> &F_i)
{
    double h=(x_opt1-0)/(N-2);

    for(int i=0; i<N-1;i++)
    {
        x_i.push_back(0+i*h);
        double F_tmp=F_2/(cosh(beta*x_i.value(i))*cosh(beta*x_i.value(i)));
        F_i.push_back(F_tmp);
    }
    x_i.push_back(l);
    F_i.push_back(F_1);
}

QVector<double> Task::f_k(double x_1, double x_2)
{
    QVector<double> f_res;
    double tmp=0.0;

    double tmp3=M*beta*a*a/(E*F_1);
    double tmp_t1=beta*x_1+M_PI/4.0;
    double tmp_t2=beta*x_2-beta*l-qAtan(tmp3)-M_PI/4.0;
    double tmp_t21=-beta*x_2+beta*l+qAtan(tmp3);
    double tmp_3=(qTan(tmp_t2))/qTan(tmp_t1);

    if(tmp_3>0) tmp=2*beta*(x_2-x_1)-qLn((qTan(tmp_t2))/qTan(tmp_t1));                          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else tmp=0.000001;

    double t=beta*l+atan(M*beta*a*a/(E*F_1));                                                   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tmp=beta*(x_2-x_1)+atanh(tan(beta*x_1))-atanh(1/(tan(t-beta*x_2)));

    f_res.push_back(tmp);

    tmp=F_1*(1-qTan(beta*x_1)*qTan(beta*x_1))-F_2*(1-1/(qTan(tmp_t21)*qTan(tmp_t21)));
    f_res.push_back(tmp);

    return f_res;
}

QVector<double> Task::matrixJacobi(double x_1, double x_2)
{
    QVector<double> J_res;

    double tmp=0.0;
    double tmp3=M*beta*a*a/(E*F_1);
    double tmp_t1=beta*x_1+M_PI/4.0;
    double tmp_t2=beta*x_2-beta*l-qAtan(tmp3)-M_PI/4.0;
    double tmp_t21=-beta*x_2+beta*l+qAtan(tmp3);

    // tmp=-2*beta+(beta)/(cos(tmp_t1)*sin(tmp_t1));
    // J_res.push_back(tmp);

    // tmp=2*beta-beta/(cos(tmp_t2)*sin(tmp_t2));
    // J_res.push_back(tmp);

    double t=beta*l+atan(M*beta*a*a/(E*F_1));

    tmp=beta*(-1+1/(cos(2*beta*x_1)));   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    J_res.push_back(tmp);                           // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    tmp=beta*(1+1/(cos(2*(t-beta*x_2))));      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    J_res.push_back(tmp);                           // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    tmp=-2*F_1*beta*sin(beta*x_1)/(cos(beta*x_1)*cos(beta*x_1)*cos(beta*x_1));
    J_res.push_back(tmp);

    tmp=2*F_1*beta*cos(tmp_t21)/(sin(tmp_t21)*sin(tmp_t21)*sin(tmp_t21));
    J_res.push_back(tmp);

    return J_res;
}

QVector<double> Task::inverseMatrixJacobi(QVector<double> J)
{
    QVector<double> J_1;
    double discriminantJ=J.value(0)*J.value(3)-J.value(1)*J.value(2);

    J_1.push_back(J.value(3)/discriminantJ);
    J_1.push_back(-J.value(1)/discriminantJ);
    J_1.push_back(-J.value(2)/discriminantJ);
    J_1.push_back(J.value(0)/discriminantJ);

    return J_1;
}

QVector<double> Task::matrixVectorMultiplication(QVector<double> matrixA, int colsA, int rowsA, QVector<double> matrixB, int colsB, int rowsB)
{
    if(colsA!=rowsB) qDebug()<<"error";
    int colsC=colsB;
    int rowsC=rowsA;
    QVector<double> matrixC;

    matrixC.push_back(matrixA.value(0)*matrixB.value(0)+matrixA.value(1)*matrixB.value(1));
    matrixC.push_back(matrixA.value(2)*matrixB.value(0)+matrixA.value(3)*matrixB.value(1));

    return matrixC;
}

double Task::findingAlpha(double x_1)
{
    double alpha = atanh(tan(beta*x_1))-beta*x_1;
    return alpha;
}
double Task::findingAlpha2(double x_2)
{
    double alpha = 0.0;
    double t=beta*l+atan(M*beta*a*a/(E*F_1));
    alpha=atanh(1.0/(tan(t-beta*x_2)))-beta*x_2;
    return alpha;
}

QVector<double> Task::decisionThreeParts()
{
    for(k_curr=k_min; k_curr<=k_max; k_curr++)
    {
        int n=searching_n_ForTheInterval_x_2(k_curr);
        x_1_curr_interval = searching_Interval_x_1(k_curr);
        x_2_curr_interval = searching_Interval_x_2(n);

        QVector<double> x_k;
        QVector<double> x_curr={0.0, 0.0};
        QVector<double> x_prev={0.0, 0.0};

        double corner =-15.0/16.0;


        while(corner<0 && x_curr.value(0)==0)
        {
            x_k.clear();
            x_curr.clear();

            x_k=findingInitialApproximationFromTheCorner(corner, true);
            if((x_k.value(0)!=0 && x_k.value(1)!=0)) x_curr = methodNewton(x_k);
            else
            {
                qDebug()<<"Метод Ньютона не нашел нужное начальное приближение x_k1 = "<<x_k<<" для интервала \nx_1="<<x_1_curr_interval<<"\nx_2="<<x_2_curr_interval<<"и угла corner = "<<corner;

                x_curr.push_back(0.0);
                x_curr.push_back(0.0);
            }
            qDebug()<<"\n\ncorner = "<<corner;
            qDebug()<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nx_curr = "<<x_curr;
            qDebug()<<"x_1="<<x_1_curr_interval<<"\nx_2="<<x_2_curr_interval;

            if(qAbs(x_curr.value(0)-x_1_curr_interval.value(0))<=10*eps || qAbs(x_curr.value(1)-x_2_curr_interval.value(0))<=10*eps || qAbs(x_curr.value(0)-x_1_curr_interval.value(1))<=10*eps || qAbs(x_curr.value(1)-x_2_curr_interval.value(1))<=10*eps)
            {
                qDebug()<<"\n\n\nКорень близок к асимптоте!!!\n\n";
                x_curr.clear();
                x_curr={0.0, 0.0};
            }
            double x_1_center=(x_1_curr_interval.value(1)-x_1_curr_interval.value(0))/2.0+x_1_curr_interval.value(0);
            if(x_curr.value(0)>x_1_center)
            {
                x_curr.clear();
                x_curr={0.0, 0.0};
            }

            corner+=1.0/16.0;
        }

        //if(x_curr.value(0)>0)          //Если нашли левую точку пересечения, запоминаем ее
        x_prev={x_curr.value(0), x_curr.value(1)};
        if(0<x_curr.value(0) &&x_curr.value(0)<l && x_1_arr.value(x_1_arr.length()-1)!=x_curr.value(0) && 0<x_curr.value(1) && x_curr.value(1)<l && x_2_arr.value(x_1_arr.length()-1)!=x_curr.value(1) && x_curr.value(0)<x_curr.value(1) && findingAlpha(x_curr.value(0))>0)
        {
            // x_1_arr.push_back(x_curr.value(0));
            // x_2_arr.push_back(x_curr.value(1));
            if(checkingTheResultsWithZeroMassAndTheTaskFromTheBook(roundingNumber(x_curr.value(0),3), roundingNumber(x_curr.value(1),3), findingAlpha(roundingNumber(x_curr.value(0),3))))
            {
                x_1_arr.push_back(roundingNumber(x_curr.value(0),3));
                x_2_arr.push_back(roundingNumber(x_curr.value(1),3));
            }
        }

        x_curr.clear();
        x_curr.push_back(0.0);
        x_curr.push_back(0.0);
        corner=-16.0;

        while(corner<0 && x_curr.value(0)==0)
        {
            x_k.clear();
            x_curr.clear();

            x_k=findingInitialApproximationFromTheCorner(corner, false);
            if((x_k.value(0)!=0 && x_k.value(1)!=0)) x_curr = methodNewton(x_k);
            else
            {
                qDebug()<<"Метод Ньютона не нашел нужное начальное приближение x_k2 = "<<x_k<<" для интервала \nx_1="<<x_1_curr_interval<<"\nx_2="<<x_2_curr_interval<<"и угла corner = "<<corner;
                x_curr.push_back(0.0);
                x_curr.push_back(0.0);
            }
            qDebug()<<"x_curr = "<<x_curr;
            qDebug()<<"corner = "<<corner;

            qDebug()<<"x_1="<<x_1_curr_interval<<"\nx_2="<<x_2_curr_interval;

            if(qAbs(x_curr.value(0)-x_1_curr_interval.value(0))<=10*eps || qAbs(x_curr.value(1)-x_2_curr_interval.value(0))<=10*eps || qAbs(x_curr.value(0)-x_1_curr_interval.value(1))<=10*eps || qAbs(x_curr.value(1)-x_2_curr_interval.value(1))<=10*eps)
            {
                qDebug()<<"\n\n\nКорень близок к асимптоте!!!\n\n";
                x_curr.clear();
                x_curr={0.0, 0.0};
            }
            else
            {
                // if(qAbs(x_curr.value(0)-x_prev.value(0))<=10*eps || qAbs(x_curr.value(1)-x_prev.value(1))<=10*eps)
                if(qAbs(x_curr.value(0)-x_prev.value(0))<=10*eps)
                {
                    qDebug()<<"\n\n\nКорень повторяется с предыдущим!!!\n\n";
                    x_curr.clear();
                    x_curr={0.0, 0.0};
                }
            }

            corner+=1.0;
        }
        x_prev.clear();
        if(0<x_curr.value(0) &&x_curr.value(0)<l && x_1_arr.value(x_1_arr.length()-1)!=x_curr.value(0) && 0<x_curr.value(1) && x_curr.value(1)<l && x_2_arr.value(x_1_arr.length()-1)!=x_curr.value(1) && x_curr.value(0)<x_curr.value(1) )
        {
            // x_1_arr.push_back(x_curr.value(0));
            // x_2_arr.push_back(x_curr.value(1));
            if(checkingTheResultsWithZeroMassAndTheTaskFromTheBook(roundingNumber(x_curr.value(0),3), roundingNumber(x_curr.value(1),3), findingAlpha(roundingNumber(x_curr.value(0),3))))
            {
                x_1_arr.push_back(roundingNumber(x_curr.value(0),3));
                x_2_arr.push_back(roundingNumber(x_curr.value(1),3));
            }

            // checkingTheResultsWithZeroMassAndTheTaskFromTheBook(x_1_arr.value(0), x_2_arr.value(0), findingAlpha(x_1_arr.value(0)));
        }
        qDebug()<<"x_curr.value(0) = "<<x_curr.value(0)<<"\nx_curr.value(1) = "<<x_curr.value(1)<<"\nalpha = "<<findingAlpha(x_curr.value(0))<<"\nalpha2 = "<<findingAlpha2(x_curr.value(1));
        qDebug()<<"\n\nx_1_curr_interval = "<<x_1_curr_interval;
        qDebug()<<"x_2_curr_interval = "<<x_2_curr_interval;

        qDebug()<<"Найденные корни";
        for(int i=0; i<x_1_arr.length(); i++)
        {
            qDebug()<<"("<<x_1_arr.value(i)<<","<< x_2_arr.value(i)<<")";
        }
    }
    qDebug()<<"x_1_arr.length() = "<<x_1_arr.length();
    if(x_1_arr.length()==0)
    {
        is_task_success=1;
        qDebug()<<"is_task_success = "<<is_task_success;
    }

    return {0.0, 0.0};

}

QVector<double> Task::methodNewton(QVector<double> x_k)
{
    QVector<double> J;
    QVector<double> J_1;
    QVector<double> f_k1;
    QVector<double> x_k2;

    QVector<double> x_k_next;
    QVector<double> x_k_prev;
    QVector<double> null_vector={0,0};
    x_k_next.push_back(x_k.value(0));
    x_k_next.push_back(x_k.value(1));

    int count=0;
    int countr=0;

    do
    {
        if(!(x_1_curr_interval.value(0)<x_k.value(0) && x_k.value(0)<x_1_curr_interval.value(1) && x_2_curr_interval.value(0)<x_k.value(1) && x_k.value(1)<x_2_curr_interval.value(1)))
        {
            qDebug()<<"\n\n!!!!Выход за пределы интервала!\n\n";
            return {0.0, 0.0};
        }
        if(qIsNaN(x_k_next.value(0)) || qIsNaN(x_k_next.value(1)))
        {
            qDebug()<<"qIsNaN(x_k_next.value(0))";
            return {0.0, 0.0};
        }
        x_k_prev.clear();

        x_k_prev.push_back(x_k.value(0));
        x_k_prev.push_back(x_k.value(1));

        x_k.clear();

        x_k.push_back(x_k_next.value(0));
        x_k.push_back(x_k_next.value(1));

        J.clear();
        J_1.clear();
        f_k1.clear();
        x_k2.clear();
        x_k_next.clear();

        J=matrixJacobi(x_k.value(0), x_k.value(1));
        J_1=inverseMatrixJacobi(J);

        f_k1=f_k(x_k.value(0), x_k.value(1));
        x_k2=matrixVectorMultiplication(J_1, 2, 2, f_k1, 1, 2);


        x_k_next.push_back(x_k.value(0)-x_k2.value(0));
        x_k_next.push_back(x_k.value(1)-x_k2.value(1));

        qDebug()<<"\n\n\nx_k = "<<x_k;
        qDebug()<<"f_k = "<<f_k1;
        qDebug()<<"J = "<<J;
        qDebug()<<"J_1 = "<<J_1;
        qDebug()<<"x_k2 = "<<x_k2;
        qDebug()<<"x_k_next = "<<QString("(%1").arg(x_k_next.value(0), 0, 'f', 8)<<QString(", %1)").arg(x_k_next.value(1), 0, 'f', 13);
        qDebug()<<"x_k_next = "<<QString("(%1, %2)").arg(x_k_next.value(0), 0, 'f', 8).arg(x_k_next.value(1), 0, 'f', 13);

        qDebug()<<"x_k_next/x_k_next; = "<<x_k_next.value(0)/x_k_next.value(0);
        count++;

        if((Norm(x_k_next, x_k)>Norm(x_k, x_k_prev) || Norm(f_k(x_k_next.value(0), x_k_next.value(1)), null_vector)>Norm(f_k(x_k.value(0), x_k.value(1)), null_vector)) && count>3)
        {
            countr++;
            if(countr>0)
            {
                qDebug()<<"Norm(x_k_next, x_k) = "<<Norm(x_k_next, x_k)<<"\nNorm(x_k, x_k_prev) = "<<Norm(x_k, x_k_prev);
                qDebug()<<"Norm(f_k(x_k_next.value(0), x_k_next.value(1)), null_vector) = "<<Norm(f_k(x_k_next.value(0), x_k_next.value(1)), null_vector)<<"\nNorm(f_k(x_k.value(0), x_k.value(1)), null_vector) = "<<Norm(f_k(x_k.value(0), x_k.value(1)), null_vector);
                qDebug()<<"Can't find roots!";
                return null_vector;
            }
        }
        if(!(x_1_curr_interval.value(0)<x_k_next.value(0) && x_k_next.value(0)<x_1_curr_interval.value(1) && x_2_curr_interval.value(0)<x_k_next.value(1) && x_k_next.value(1)<x_2_curr_interval.value(1)))
        {
            qDebug()<<"\n\n!!!!Выход за пределы интервала222!\n\n";
            return {0.0, 0.0};
        }

        if(qIsNaN(x_k_next.value(0)) || qIsNaN(x_k_next.value(1)))
        {
            qDebug()<<"qIsNaN(x_k_next.value(0))2";
            return {0.0, 0.0};
        }
    }while(Norm(x_k_next, x_k)>eps && count<3000);
    qDebug()<<"count = "<<count;
    qDebug()<<"f = "<<f_k(x_k_next.value(0), x_k_next.value(1));
    qDebug()<<"Norm(x_k_next, x_k) = "<<Norm(x_k_next, x_k)<<"\neps = "<<eps;

    qDebug()<<"Numbers of iterations = "<<count;

    // if(x_1_interval.value(0)<x_k_next.value(0) && x_k_next.value(0)<x_1_interval.value(1) && x_2_interval.value(0)<x_k_next.value(1) && x_k_next.value(1)<x_2_interval.value(1) && 0<x_k_next.value(0) && 0<x_k_next.value(1)) return {x_1_curr, x_2_curr};

    // return {0.0, 0.0};
    return x_k_next;
}




void Task::fillingInVectors3(QVector<double> &x_i, QVector<double> &F_i, QVector<double> x_opt)
{
    double h=(x_opt.value(1)-x_opt.value(0))/(N-3);

    x_i.push_back(0);
    F_i.push_back(F_2);

    double alpha=0.0;
    alpha=atanh(tan(beta*x_opt.value(0)))-beta*x_opt.value(0);
    qDebug()<<"alpha2 = "<<alpha;
    qDebug()<<"F_2*(cosh(beta*x_opt.value(0)+alpha)*cosh(beta*x_opt.value(0)+alpha)) = "<<F_2*(cosh(beta*x_opt.value(0)+alpha)*cosh(beta*x_opt.value(0)+alpha));

    x_i.push_back(x_opt.value(0)-h);
    F_i.push_back(F_2);

    for(int i=0; i<N-2;i++)
    {
        x_i.push_back(x_opt.value(0)+i*h);

        // qDebug()<<"cosh(beta*x_opt.value(1)+alpha)*cosh(beta*x_opt.value(1)+alpha) = "<<cosh(beta*x_opt.value(1)+alpha)*cosh(beta*x_opt.value(1)+alpha);
        // qDebug()<<"cosh(beta*x_i.value(i)+alpha)*cosh(beta*x_i.value(i)+alpha) = "<<cosh(beta*x_i.value(i)+alpha)*cosh(beta*x_i.value(i)+alpha);

        // double F_tmp=F_1*(cosh(beta*x_opt.value(1)+alpha)*cosh(beta*x_opt.value(1)+alpha))/(cosh(beta*x_i.value(i)+alpha)*cosh(beta*x_i.value(i)+alpha));
        double F_tmp=F_2*(cosh(beta*x_opt.value(0)+alpha)*cosh(beta*x_opt.value(0)+alpha))/(cosh(beta*x_i.value(i+2)+alpha)*cosh(beta*x_i.value(i+2)+alpha));
        F_i.push_back(F_tmp);

    }
    x_i.push_back(l);
    F_i.push_back(F_1);

    qDebug()<<"x_i = "<<x_i;
    qDebug()<<"F_i = "<<F_i;
    qDebug()<<"l = "<<l;
    qDebug()<<"x_opt = "<<x_opt;
}

int Task::searching_n_ForTheInterval_x_2(int k)
{
    // qDebug()<<"floor(-1.234) = "<<floor(-1.234);
    // qDebug()<<"floor(1.234) = "<<floor(1.234);
    // qDebug()<<"ceil(-1.234) = "<<ceil(-1.234);
    // qDebug()<<"ceil(1.234) = "<<ceil(1.234);

    double t=beta*l+atan(M*beta*a*a/(E*F_1));

    int n=ceil(1/2 + k -t/M_PI);
    return n;
}

QVector<double> Task::searching_Interval_x_1(int k)
{
    QVector<double> x_1_interval;
    double tmp = (3*M_PI/4.0 +M_PI*k)/beta;

    x_1_interval.push_back(tmp);
    tmp = (M_PI/4.0 +M_PI*(k+1))/beta;
    x_1_interval.push_back(tmp);
    return x_1_interval;
}

QVector<double> Task::searching_Interval_x_2(int n)
{
    double t=beta*l+atan(M*beta*a*a/(E*F_1));
    QVector<double> x_2_interval;

    double tmp = (t +M_PI/4.0 +M_PI*n)/beta;
    x_2_interval.push_back(tmp);

    tmp=(t + 3*M_PI/4.0 +M_PI*n)/beta;
    x_2_interval.push_back(tmp);

    return x_2_interval;
}




QVector<double> Task::findingInitialApproximationFromTheCorner(double corner, bool is_left_root)
{
    int count=0;
    // k_curr=-1;
    int k=k_curr;
    // k=0;
    int n=searching_n_ForTheInterval_x_2(k);

    // QVector<double> x_1_interval = searching_Interval_x_1(k);
    // QVector<double> x_2_interval = searching_Interval_x_2(n);

    QVector<double> x_1_interval = x_1_curr_interval;
    QVector<double> x_2_interval = x_2_curr_interval;
    double t=beta*l+atan(M*beta*a*a/(E*F_1));

    // double s=-1.0/16.0;
    double s=corner;

    double s1=0.0;
    double x_1_prev;
    double x_1_curr;
    double x_1_next;

    if(is_left_root)
    {
        s1=-s*(x_1_interval.value(0)+eps)+(x_2_interval.value(1)-eps);
        x_1_curr=x_1_interval.value(0)+eps;
    }
    else
    {
        s1=-s*(x_1_interval.value(1)-eps)+(x_2_interval.value(0)+eps);
        x_1_curr=x_1_interval.value(1)-eps;

        s=-(17.0+corner)/16.0;
        double x1_tmp=(x_1_interval.value(1)-x_1_interval.value(0))/2+x_1_interval.value(0);
        s1=-s*(x1_tmp+eps)+(x_2_interval.value(1)-eps);
        x_1_curr=x1_tmp+eps;
        // x_1_curr=x_1_interval.value(1)-eps;
        // x_1_curr=-(x_1_interval.value(1)-x_1_interval.value(0))/2+x_1_interval.value(1);
        // if(s1>0) x_1_curr=x_1_interval.value(1)-eps;
        // else x_1_curr=x_1_curr=(x_1_interval.value(1) -x_1_interval.value(0))/2.0+x_1_interval.value(0)+eps;
        // x_1_curr=x_1_interval.value(0)+eps;
        // x_1_curr=(x_1_interval.value(1) -x_1_interval.value(0))/2.0+x_1_interval.value(0)-eps;
    }

    double t1=t-beta*s1;
    x_1_next=x_1_curr;
    double f1, f1_prev;
    double df1dx1;
    int countr=0;
    int countr_count=0;


    do
    {
        qDebug()<<"x_1_curr = "<<x_1_curr;
        x_1_prev=x_1_curr;
        x_1_curr=x_1_next;
        df1dx1=beta*(s-1+1.0/(cos(2*beta*x_1_curr))+s/(cos(2*(t1-beta*s*x_1_curr))));
        f1=beta*(s*x_1_curr +s1-x_1_curr)+atanh(tan(beta*x_1_curr))-atanh(1.0/(tan(t-beta*s*x_1_curr-beta*s1)));
        f1_prev=f1;
        x_1_next=x_1_curr-f1/df1dx1;

        count++;

        // if((qAbs(x_1_next-x_1_curr)>qAbs(x_1_curr-x_1_prev) || qAbs(f1)>qAbs(f1_prev)) && count>3)
        if((qAbs(f1)>qAbs(f1_prev)) && count>3)
        {
            if(countr_count==count-1)
            {
                qDebug()<<"count = "<<count;
                countr++;
            }
            else countr=0;

            if(countr>3)
            {
                qDebug()<<"Can't find roots with corner!";
                return {0.0, 0.0};
            }
            countr_count=count;
        }
    }while(qAbs(x_1_next-x_1_curr)>=eps && count<55);

    qDebug()<<"initial_x_1 = "<<x_1_curr;
    double x_2_curr=s*x_1_curr+s1;
    qDebug()<<"("<<x_1_curr<<","<<x_2_curr<<")";

    if(x_1_interval.value(0)<x_1_curr && x_1_curr<x_1_interval.value(1) && x_2_interval.value(0)<x_2_curr && x_2_curr<x_2_interval.value(1)) return {x_1_curr, x_2_curr};
    return {0.0, 0.0};
}

double Task::Norm(QVector<double> x_i, QVector<double> x_i_prev)
{
    double res=0.0;
    int length=x_i.length();

    for(int i=0; i<length; i++)
    {
        res+=(x_i.value(i)-x_i_prev.value(i))*(x_i.value(i)-x_i_prev.value(i));
    }
    res=qSqrt(res);
    return res;
}

bool Task::checkingTheResultsWithZeroMassAndTheTaskFromTheBook(double x_1, double x_2, double alpha)
{
    // double equation1 = tanh(beta*x_2+alpha)-1.0/(tan(beta*(l-x_2)));
    double s= beta*l+atan(M*beta*a*a/(E*F_1));
    double equation1 = tanh(beta*x_2+alpha)-1.0/(tan(s-beta*(x_2)));
    // double equation2 = 1/tan(s-beta*x_2)-tanh(beta*x_1+alpha);
    double equation2 = tan(beta*x_1)-tanh(beta*x_1+alpha);
    double equation3=F_1*cosh(beta*x_2+alpha)*cosh(beta*x_2+alpha)-F_2*cosh(beta*x_1+alpha)*cosh(beta*x_1+alpha);

    QVector<double> J=matrixJacobi(x_1, x_2);

    qDebug()<<QString("equation1 = %1\nequation1 = %2\nequation1 = %3\n").arg(equation1, 0, 'f', 16).arg(equation2, 0, 'f', 16).arg(equation3, 0, 'f', 16);
    // qDebug()<<QString("equation1 = %1\nequation1 = %2\nequation1 = %3\n").arg(10*eps, 0, 'f', 16).arg(10*eps, 0, 'f', 16).arg((J.value(2)+J.value(3))*eps, 0, 'f', 16);
    // if(fabs(equation1)<=delta && fabs(equation2)<=delta && fabs(equation3)<=(J.value(2)+J.value(3))*eps)
    if(fabs(equation1)<=delta && fabs(equation2)<=delta && fabs(equation3)<=delta)
    {
        qDebug()<<"SUCCESS!";
        return true;
    }
    else qDebug()<<"NOOO";
    return false;
}


