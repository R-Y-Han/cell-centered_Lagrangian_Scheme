#ifndef INITIAL_H
#define INITIAL_H


//��������һЩ��ֵ�Ͳ������˴�x��ʾ������y��ʾ�Ƕ�
static const double pi = 3.1415926;

static const int nx=5;  //x���򻮷ָ���
static const int ny=4;  //y���򻮷ָ���
static const int n= 2 * nx * (ny+1) + 1; //�ܵĽڵ��������Ϊ��0��ʼ��

static const double r_1 = 0.1, r_2 = 0.11;
static const double x = r_2, y=pi / 2;    //x���򳤶ȣ�y���򳤶�

static const double hx_1= (double) r_1/nx, hx_2 = (double) (r_2 - r_1)/nx,  hy = (double) y/ny;  //x���򲽳���y���򲽳�

static const double rho_1 = 0.01, rho_2 = 1, P_1 = 2.5e9, P_2 = 2.5e11;  //��ʼ״̬�ܶȡ�ѹǿ��ֵ
static const double gamma_0 = (double) 5/3;    //״̬���̲���
static const double vn=0;    //�߽��Ϸ����ٶ�Ϊ0

static const double Ce = 0.3, Cv = 0.1, Cm = 1.01;   //CFL������ϵ��

static const double a_0 = 0;

#endif