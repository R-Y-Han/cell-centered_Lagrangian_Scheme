#ifndef INITIAL_H
#define INITIAL_H


//用来声明一些初值和步长，此处x表示极径，y表示角度
static const double pi = 3.1415926;

static const int nx=10;  //x方向划分个数
static const int ny=40;  //y方向划分个数
static const int n= (nx+1) * (ny+1); //总的节点个数，因为从0开始数

static const double r_1 = 0.9, r_2 = 1; //内外边界
static const double x = r_2 - r_1, y=pi ;    //x方向长度，y方向长度

static const double hx= (double) x/nx,  hy = (double) y/ny;  //x方向步长，y方向步长

static const double rho_1 = 1e-3, rho_2 = 1e-2, P_1 = 0.1, P_2  = 10;  //初始状态密度、压强的值
static const double s = 1e5;    //熵
static const double t_end = sqrt(0.5 * (r_2*r_2 - r_1 * r_1)/(2*s*(rho_2 - rho_1))); //终止时间
static const double gamma_0 = (double) 2;    //状态方程参数
static const double vn=0;    //边界上法向速度为0

static const double Ce = 0.1, Cv = 0.1, Cm = 1.01;   //CFL的三个系数

//扰动相关系数
static const int mode = 8;
static const double xi_1 = 0, xi_2 = 1;
static const double a_0 = 1e-3;
static double A_1, A_2;

#endif