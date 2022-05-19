#ifndef INITIAL_H
#define INITIAL_H


//��������һЩ��ֵ�Ͳ������˴�x��ʾ������y��ʾ�Ƕ�
static const double pi = 3.1415926;

static const int nx=10;  //x���򻮷ָ���
static const int ny=40;  //y���򻮷ָ���
static const int n= (nx+1) * (ny+1); //�ܵĽڵ��������Ϊ��0��ʼ��

static const double r_1 = 0.9, r_2 = 1; //����߽�
static const double x = r_2 - r_1, y=pi ;    //x���򳤶ȣ�y���򳤶�

static const double hx= (double) x/nx,  hy = (double) y/ny;  //x���򲽳���y���򲽳�

static const double rho_1 = 1e-3, rho_2 = 1e-2, P_1 = 0.1, P_2  = 10;  //��ʼ״̬�ܶȡ�ѹǿ��ֵ
static const double s = 1e5;    //��
static const double t_end = sqrt(0.5 * (r_2*r_2 - r_1 * r_1)/(2*s*(rho_2 - rho_1))); //��ֹʱ��
static const double gamma_0 = (double) 2;    //״̬���̲���
static const double vn=0;    //�߽��Ϸ����ٶ�Ϊ0

static const double Ce = 0.1, Cv = 0.1, Cm = 1.01;   //CFL������ϵ��

//�Ŷ����ϵ��
static const int mode = 8;
static const double xi_1 = 0, xi_2 = 1;
static const double a_0 = 1e-3;
static double A_1, A_2;

#endif