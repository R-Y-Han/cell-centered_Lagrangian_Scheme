#ifndef INITIAL_H
#define INITIAL_H


//用来声明一些初值和步长，此处x表示极径，y表示角度
static const double pi = 3.1415926;

static const int nx=5;  //x方向划分个数
static const int ny=4;  //y方向划分个数
static const int n= 2 * nx * (ny+1) + 1; //总的节点个数，因为从0开始数

static const double r_1 = 0.1, r_2 = 0.11;
static const double x = r_2, y=pi / 2;    //x方向长度，y方向长度

static const double hx_1= (double) r_1/nx, hx_2 = (double) (r_2 - r_1)/nx,  hy = (double) y/ny;  //x方向步长，y方向步长

static const double rho_1 = 0.01, rho_2 = 1, P_1 = 2.5e9, P_2 = 2.5e11;  //初始状态密度、压强的值
static const double gamma_0 = (double) 5/3;    //状态方程参数
static const double vn=0;    //边界上法向速度为0

static const double Ce = 0.3, Cv = 0.1, Cm = 1.01;   //CFL的三个系数

static const double a_0 = 0;

#endif