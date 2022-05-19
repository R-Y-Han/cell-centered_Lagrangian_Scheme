#ifndef INITIAL_PO_H
#define INITIAL_PO_H


//用来声明一些初值和步长

static const int nx=50;  //x方向划分个数
static const int ny=25;  //y方向划分个数
static const int n= nx*(ny+1)+1; //总的节点个数，因为从0开始数，所以加一
static const double pi = 3.1415926;

static const double x = 1, y=pi/2;    //x方向长度，y方向长度，在x/2处状态有间断

static const double hx= (double) x/nx,  hy = (double) y/ny;  //x方向步长，y方向步长

static const double rho_0 = 1, P_0 = 1e-9, u_0 = 0;  //初始状态活密度、压强、速度 的值
static const double gamma_0 = (double) 5/3;    //状态方程参数
static const double vn =0;    //边界上法向速度为0

static const double Ce = 0.3, Cv = 0.1, Cm = 1.01;   //CFL的三个系数


#endif