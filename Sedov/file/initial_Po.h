#ifndef INITIAL_PO_H
#define INITIAL_PO_H


//��������һЩ��ֵ�Ͳ���

static const int nx=50;  //x���򻮷ָ���
static const int ny=25;  //y���򻮷ָ���
static const int n= nx*(ny+1)+1; //�ܵĽڵ��������Ϊ��0��ʼ�������Լ�һ
static const double pi = 3.1415926;

static const double x = 1, y=pi/2;    //x���򳤶ȣ�y���򳤶ȣ���x/2��״̬�м��

static const double hx= (double) x/nx,  hy = (double) y/ny;  //x���򲽳���y���򲽳�

static const double rho_0 = 1, P_0 = 1e-9, u_0 = 0;  //��ʼ״̬���ܶȡ�ѹǿ���ٶ� ��ֵ
static const double gamma_0 = (double) 5/3;    //״̬���̲���
static const double vn =0;    //�߽��Ϸ����ٶ�Ϊ0

static const double Ce = 0.3, Cv = 0.1, Cm = 1.01;   //CFL������ϵ��


#endif