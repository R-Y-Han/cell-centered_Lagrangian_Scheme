#ifndef INITIAL_H
#define INITIAL_H


//用来声明一些初值和步长

static const int nx=100;  //x方向划分个数
static const int ny=10;  //y方向划分个数
static const int n= (nx+1)*(ny+1); //总的节点个数，因为从0开始数，所以加一

static const double x = 1, y=0.1;    //x方向长度，y方向长度，在x/2处状态有间断

static const double hx= (double) x/nx,  hy = (double) y/ny;  //x方向步长，y方向步长

static const double rho_l = 1, P_l = 1, u_l = 0;  //初始状态活塞左端 密度、压强、速度 的值
static const double rho_r = 0.125, P_r = 0.1, u_r = 0;   //处置状态活塞右端的值
static const double gamma_l = (double) 7/5, gamma_r = (double) 5/3;    //状态方程参数
static const double vn=0;    //边界上法向速度为0

static const double Ce = 0.3, Cv = 0.1, Cm = 1.01;   //CFL的三个系数


#endif