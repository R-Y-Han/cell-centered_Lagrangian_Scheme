#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "initial.h"

using namespace std;

//************下面定义需要用到的数组***********//

//关于节点的量

static double loc[2][n+1][2]; //loc[a][qm][k]为第(i,j)点在第ta时刻的第k个分量，ta=tn或tn+1
static double nodev[n+1][2];    //nodev[qm][k]为第(i,j)个点的第k个速度分量，根据算法只需要一个数组存储
static int boundary[n+1];   //等于0代表是内点，1是边上边界点，2是弧上边界点

//关于网格的量

static double tau[2][nx][ny];   //tau[a][h][l]代表第(h,l)个网格的比容，ta时刻
static double cellv[2][nx][ny][2];  //cellv[a][h][l][k]为第(h,l)个网格的第k个速度分量，ta时刻
static double E[2][nx][ny]; //E[a][h][l]为(h,l)个网格的总能量，ta时刻
static double cellP[2][nx][ny]; //cellP[a][h][l]为(h,l)网格的压力，ta时刻
static double m[nx][ny];    //第(h,l)个网格的质量，守恒
static double e[2][nx][ny]; //第(h,l)个网格的内能，ta时刻
static double c[2][nx][ny]; //第(h,l)个网格的声速，ta时刻
static double A[2][nx][ny]; //第(h,l)个网格的面积，ta时刻
static double gamma[nx][ny];    //第(h,l)个网格的状态方程参数，守恒

//关于划分的量

vector<int> neighbornode[n+1];   //第(i,j)个节点的所有相邻节点
vector<int> neighborcell[n+1];   //第(i,j)个网格的所有邻域
vector<double> L[2][n+1];    //第(i,j)个节点的第k条边的边长
vector<double> Nx[2][n+1];   //第(i,j)个节点的第k条边外方向x分量
vector<double> Ny[2][n+1];   //第(i,j)个节点的第k条边外方向y分量
vector<double> Pstar[n+1][2];    //第(i,j)个节点第k条边两边的压力
//******************数组定义完毕**************//

enum time {tlast, tnext};   //定义枚举类型，tlast为上一个时间点，tnest为下一个时间点，分别代表0，1

//**************需要的函数声明****************//
void initialmesh(); //初始化

void internalsolver(int t_l, int qm);  //计算内点的速度
void getinternalPstar(int t_l, int qm);   //计算内点边上的压力
void boundarysolverNP(int t_l, int qm);  //计算边界点（不平行）的速度
void boundarysolverP(int t_l, int qm);    //计算边界点（平行）的速度
void anglesolver1(int t_l, int qm); //(nx,0)点的速度
void anglesolver2(int t_l, int qm); //(nx,ny)点的速度
void boundarywithP(int t_l, int qm);    //边界条件由压力给出的速度
void boundarysolver2(int t_l, int qm);  //计算边界点类型为2时的速度
void getboundaryPstar(int t_l, int qm);   //计算边界点上边上的压力
void nodalsolver(int t_l); //计算上一时刻网格节点的速度，根据节点性质选择上面三个函数中的一个
double choosedt(int t_l, double dtn);   //计算时间步长
void updategeometric(int t_l, int t_n, double dt);  //更新新时间点的几何量
void updatephysicalvariables(int t_l, int t_n, double dt);   //更新新时间点的物理量
void getstates(int t_n);    //根据状态方程求解剩下变量

//***************函数声明完毕****************//

int main()
{
    
    int i, j, h, l, q;
    time tl, tn, ttemp;
    double t = 0;
    double dt = 1000;
    double temp;
    
    initialmesh();

    tl = tlast;
    tn = tnext;
    /*
    for (q=0; q<=n; q++)
    {
        if (q == 0)
        {
            i = 0;
            j = 0;
        }
        else{
            i = (q-1) / (ny + 1) + 1;
            j = (q-1) % (ny + 1);
        }

        cout<<q<<"(";
        for (h=0; h<neighborcell[q].size(); h++)
        {
            cout<<neighborcell[q][h]<<"\t";
        }
        cout<<")"<<"\t";
        if ( q % (ny + 1) == 0 )
        {
            cout<<endl;
        }

    }*/
    /*****************************************
    {
        fstream test;
        remove("cell-centered_Lagrangian_Scheme\\Noh\\output\\mesh.plt");
        test.open("cell-centered_Lagrangian_Scheme\\Noh\\output\\mesh.plt", ios :: out | ios :: app);
        test<<"ZONE"<<"\t"<<"i="<<4<<", "<<"j="<<3<<", F = POINT"<<endl;
        for (q = 0; q<=n; q++)
        {
            if (q == 0)
            {
                i = 0;
                j = 0;
            }
            else{
                i = (q - 1) / (ny + 1) + 1;
                j = (q - 1) % (ny + 1);
            }

            cout<<i<<"\t"<<j<<"\t"<<cos(j*hy)*hx*i<<"\t"<<sin(j*hy)*hx*i<<endl;
            test<<"\t"<<loc[tl][q][0]<<"\t"<<loc[tl][q][1]<<endl;
        }
        test.close();
    }
    //******************************************/
    
    fstream fv, frho, fp, fe;   //四个流传输

    //文件路径
    const char* filev = "cell-centered_Lagrangian_Scheme\\Noh\\output\\V.txt";    //记录单元速度
    const char* filerho = "cell-centered_Lagrangian_Scheme\\Noh\\output\\rho.txt";    //记录单元密度
    const char* filep = "cell-centered_Lagrangian_Scheme\\Noh\\output\\P.txt";    //记录单元压力
    const char* filee = "cell-centered_Lagrangian_Scheme\\Noh\\output\\e.txt";    //记录单元内能

    remove(filev);
    remove(filerho);
    remove(filep);
    remove(filee);
    
    do{
        //先记录数据，结果太大了不留了（600M）

        /*t
        //记录速度
        {
            fv.open(filev,ios::out| ios::app);
            fv<<"t = "<<t<<endl;
            fv<<"\\delta t = "<<dt<<endl;
        
            for (h=0; h<nx; h++)
            {
                for (l=0; l<ny; l++)
                {
                    fv<<"("<<cellv[tl][h][l][0]<<","<<cellv[tl][h][l][1]<<")"<<"\t";
                }
                fv<<endl<<endl;
            }

            fv.close();
        }

        //记录单元密度
        {
            frho.open(filerho,ios::out| ios::app);
            frho<<"t = "<<t<<endl;
            frho<<"\\delta t = "<<dt<<endl;
        
            for (h=0; h<nx; h++)
            {
                for (l=0; l<ny; l++)
                {
                    temp = 1 / tau[tl][h][l];
                    frho<<temp<<"\t";
                }
                frho<<endl<<endl;
            }

            frho.close();
        }

        //记录单元压力
        {
            fp.open(filep,ios::out| ios::app);
            fp<<"t = "<<t<<endl;
            fp<<"\\delta t = "<<dt<<endl;
        
            for (h=0; h<nx; h++)
            {
                for (l=0; l<ny; l++)
                {
                    fp<<cellP[tl][h][l]<<"\t";
                }
                fp<<endl<<endl;
            }

            fp.close();
        }

        //记录网格内能
        {
            fe.open(filee,ios::out| ios::app);
            fe<<"t = "<<t<<endl;
            fe<<"\\delta t = "<<dt<<endl;
        
            for (h=0; h<nx; h++)
            {
                for (l=0; l<ny; l++)
                {
                    fe<<e[tl][h][l]<<"\t";
                }
                fe<<endl<<endl;
            }

            fe.close();
        }
*/
        //进行运算
        nodalsolver(tl);

        dt = choosedt(tl, dt);
        t = t+ dt;

        updategeometric(tl,tn,dt);
        updatephysicalvariables(tl,tn,dt);
        getstates(tn);

        ttemp = tl;
        tl = tn;
        tn = ttemp; //原来tl时刻的网格求解得tn时刻的网格，下一次计算时tn就为上一时刻的网格
        cout<<t<<"\t"<<dt<<endl<<endl;
    }while(t < 0.6);

    //下面开始画图
    {
        //文件名
        const char* plotv = "cell-centered_Lagrangian_Scheme\\Noh\\output\\V.plt";    //单元速度
        const char* plotrho = "cell-centered_Lagrangian_Scheme\\Noh\\output\\rho.plt";    //单元密度
        const char* plotp = "cell-centered_Lagrangian_Scheme\\Noh\\output\\P.plt";    //单元压力
        const char* plote = "cell-centered_Lagrangian_Scheme\\Noh\\output\\e.plt";    //单元内能
        const char* point = "cell-centered_Lagrangian_Scheme\\Noh\\output\\mesh.plt";

        remove(plotv);
        remove(plotrho);
        remove(plotp);
        remove(plote);
        remove(point);

        fstream fpoint;
    
        //先画速度，因为是sod shock，所以y向速度为0，只需画x向速度
        {
            fv.open(plotv, ios :: out | ios :: app);
            fv<<"VARIABLES="<<"X"<<","<<"V"<<endl;
            for (h=0; h<nx; h++)
            {
                if ( h == 0)
                {
                    q = 0;
                }
                else{
                    q = (h - 1) * (ny + 1) + 2;
                }

                fv<<"\t"<<loc[tl][q][0]<<"\t"<<cellv[tl][h][1][0]<<endl;
            }
            fv.close();
        }

        //再画密度，由于sod shock是一维问题，所以只要记录x向，下同
        {
            frho.open(plotrho, ios :: out | ios :: app);
            frho<<"VARIABLES="<<"X"<<","<<"density"<<endl;
            for (h=0; h<nx; h++)
            {
                if ( h == 0)
                {
                    q = 0;
                }
                else{
                    q = (h - 1) * (ny + 1) + 2;
                }
                temp = 1 / tau[tl][h][1];

                frho<<"\t"<<loc[tl][q][0]<<"\t"<<temp<<endl;
            }
            frho.close();
        }

        //再画压力
        {
            fp.open(plotp, ios :: out | ios :: app);
            fp<<"VARIABLES="<<"X"<<","<<"P"<<endl;
            for (h=0; h< nx; h++)
            {
                if ( h == 0)
                {
                    q = 0;
                }
                else{
                    q = (h - 1) * (ny + 1) + 2;
                }

                fp<<"\t"<<loc[tl][q][0]<<"\t"<<cellP[tl][h][1]<<endl;
            }
            fp.close();
        }

        //再画内能
        {
            fe.open(plote, ios :: out | ios :: app);
            fe<<"VARIABLES="<<"X"<<","<<"e"<<endl;
            for (h=0; h< nx; h++)
            {
                if ( h == 0)
                {
                    q = 0;
                }
                else{
                    q = (h - 1) * (ny + 1) + 2;
                }

                fe<<"\t"<<loc[tl][q][0]<<"\t"<<e[tl][h][1]<<endl;
            }
            fe.close();
        }

        //再画网格
        {
            fpoint.open(point, ios :: out | ios :: app);
            fpoint<<"VARIABLES="<<"X"<<","<<"Y"<<endl;
            fpoint<<"ZONE N="<<n+1<<", E="<<nx*ny<<", F=FEPOINT, ET=QUADRILATERAL"<<endl;
            for (q=0; q<=n; q++)
            {
                fpoint<<"\t"<<loc[tl][q][0]<<"\t"<<loc[tl][q][1]<<endl;
            }
            double ctemp[4];
                
            for (h=0; h<nx; h++)
            {
                for (l=0; l<ny; l++)
                {
                    if (h == 0)
                    {
                        ctemp[0] = 0;
                        ctemp[1] = l + 1;
                        ctemp[2] = l+2;
                        ctemp[3] = 0;
                    }
                    else{
                        ctemp[0] = (h-1) * (ny+1) + l +1;
                        ctemp[1] = h * (ny+1) + l + 1;
                        ctemp[2] = h * (ny+1) + l + 2;
                        ctemp[3] = (h-1) * (ny+1) + l +2;
                    }
                    
                    for (q=0; q<4; q++)
                    {
                        fpoint<<"\t"<<ctemp[q]+1;
                    }
                    fpoint<<endl;
                }
            }
            fpoint.close();
        }
    }
    
    system("pause");
}

//***********注意我在这里初始化网格速度的时候和文章不同***********//
void initialmesh()
{
    int i, j, h, l, q, qt;
    double temp, temp1, temp2;

    //先初始化节点相关量
    for (q=0; q<=n; q++)
    {
        //节点坐标
        if (q == 0)
        {
            i = 0;
            j = 0;
        }
        else{
            i = (q - 1) / (ny + 1) + 1;
            j = (q - 1) % (ny + 1);
        }

            loc[0][q][0] = (i * hx) * cos(j * hy);   //x坐标
            loc[0][q][1] = (i * hx) * sin(j * hy);   //y坐标

            //判断是否是边界点
            if (i == 0)
            {
                boundary[q] = 1; //顶点，此时边界条件由法向速度给出
            }
            else if ( i == nx)
            {
                boundary[q] = 2; //弧长上的点，此时bc由压力给出
            }
            else if (j == 0 || j == ny)
            {
                boundary[q] = 1; //边上的点，此时bc由法向速度给出
            }
            else{
                boundary[q] = 0;
            }   //边界点初始化完毕
    }
    //节点初始化完毕

    //初始化网格量
    for (h=0; h<nx; h++)
    {
        for (l=0; l<ny; l++)
        {
            
            {
                tau[0][h][l] = 1/ rho_0;

                if (h == 0)
                {
                    cellv[0][h][l][0] = 0;
                    cellv[0][h][l][1] = 0;
                }
                else{
                    cellv[0][h][l][0] = - cos((l+0.5) * hy);
                    cellv[0][h][l][1] = - sin((l+0.5) * hy);
                }
                

                cellP[0][h][l] = P_0;

                gamma[h][l] = gamma_0;

                e[0][h][l] = cellP[0][h][l] * tau[0][h][l] / (gamma[h][l] - 1);  //内能

                E[0][h][l] = 0.5 * (cellv[0][h][l][0] * cellv[0][h][l][0] + cellv[0][h][l][1] * cellv[0][h][l][1]);
                E[0][h][l] = E[0][h][l] + e[0][h][l];   //总能量
            
                c[0][h][l] = (gamma[h][l] - 1) * e[0][h][l];
                c[0][h][l] = sqrt(c[0][h][l] ); //声速

                A[0][h][l] = 0.5 * ( (h+1) * hx ) * ((h+1) * hx) * hy;
                A[0][h][l] = A[0][h][l] - 0.5 * ( h * hx ) * (h * hx) * hy;   //面积
                m[h][l] = A[0][h][l] / tau[0][h][l];   //质量
            }
            
        }
    }
    //网格量初始化完毕

    //初始化划分量

    //先对顶点进行初始化
    {
        //先分配邻点
        for (l = 0; l<=ny; l++)
        {
            neighbornode[0].push_back(l+1);
        }

        //再分配相邻单元
        for (l=0; l<ny; l++)
        {
            neighborcell[0].push_back(l);
        }
        neighborcell[0].push_back(-1);

        //再计算长度
        {
            for (l=0; l<neighbornode[0].size(); l++)
            {
                q = neighbornode[0][l];
                temp = (loc[0][q][0] - loc[0][0][0])*(loc[0][q][0] - loc[0][0][0]);
                temp = temp + (loc[0][q][1] - loc[0][0][1])*(loc[0][q][1] - loc[0][0][1]);
                temp = sqrt(temp);

                L[0][0].push_back(temp);
                L[1][0].push_back(0);
            }
        }

        //再计算外法向量
        {
            for (l=0; l<neighbornode[0].size(); l++)
            {
                q = neighbornode[0][l];
                temp1 = loc[0][0][1] - loc[0][q][1];
                temp2 = loc[0][q][0] - loc[0][0][0];
                temp1 = temp1 / L[0][0][l];
                temp2 = temp2 / L[0][0][l];

                Nx[0][0].push_back(temp1);
                Ny[0][0].push_back(temp2);
                Nx[1][0].push_back(0);
                Ny[1][0].push_back(0);

                Pstar[0][0].push_back(0);
                Pstar[0][1].push_back(0);
            }
        }
    }

    //再对剩下节点初始化
    int* ntemp = new int [4];
    int* ctemp = new int [4];
    for (q=1; q<=n; q++)
    {

        i = (q - 1) / (ny + 1) + 1;
        j = (q - 1) % (ny + 1);
        
        //先找每个点的邻点
        {
            if (j-1 < 0)    //(i,j-1)
            {
                ntemp[0] = -1;
            }
            else{
                ntemp[0] = (i - 1) * (ny + 1) + j - 1 + 1;
            }
            if (i+1 > nx)   //(i+1,j)
            {
                ntemp[1] = -1;
            }
            else{
                ntemp[1] = (i+1 - 1) * (ny + 1 ) + j + 1;
            }
            if (j+1 > ny)   //(i,j+1)
            {
                ntemp[2] = -1;
            }
            else{
                ntemp[2] = (i - 1) * (ny + 1) + j + 1 + 1;
            }
            if (i-1 < 0)    //(i-1,j)
            {
                ntemp[3] = -1;
            }
            else{
                if ( i == 1)
                {
                    ntemp[3] = 0;
                }
                else{
                    ntemp[3] = (i-1 - 1) * (ny + 1) + j + 1;
                }
            }
        }
        
        //再按逆时针邻点连续出现排列
        int a=0, b=0, bb;
        for (a=0; a<4; a++)
        {
            b = a+1;
            if (b>3)
            {
                b = b - 4;
            }

            if (ntemp[a] < 0 && ntemp[b] >= 0)
            {
                break;
            }
        }
        //b为邻点连续出现时第一个邻点的编号
        bb = b;

        b = bb;
        //记录第(i,j)个节点的邻点编号
        for (a=0 ; a<4; a++)
        {
            if (ntemp[b] >= 0)
            {
                neighbornode[q].push_back(ntemp[b]);
            }
                
            b++;
            if (b>3)
            {
                b = b - 4;
            }
        }
        //记录邻点编号完毕

        //记录邻域编号
        {
            ctemp[0] = i * ny + j - 1;
            ctemp[1] = i * ny + j;
            if (i == 1)
            {
                ctemp[2] = j;
                ctemp[3] = j - 1;
            }
            else{
                ctemp[2] = (i-1) * ny + j;
                ctemp[3] = (i-1) * ny + j - 1;
            }
        }
            
        //判断邻域类型
        for (a = 0; a < 4; a++)
        {
            b = a + 1; 
            if (b>3)
            {
                b = b - 4;
            }

            if (ntemp[a] < 0)
            {
                ctemp[a] = -2;  //-2表示不需要考虑
            }
            else{
                if (ntemp[b] < 0 )
                {
                    ctemp[a] = -1;  //-1表示虚网格
                }
            }
        }

        //记录邻域
        b = bb;
        for (a=0; a< 4; a++)
        {
            if (ctemp[b] > -2)
            {
                neighborcell[q].push_back(ctemp[b]);
            }
            b++;
            if (b>3)
            {
                b = b-4;
            }
        }
        //记录邻域完毕

        //计算每个点到相邻节点的长度和法向量
        for (bb = 0; bb < neighbornode[q].size(); bb++)
        {
            qt = neighbornode[q][bb];

            temp1 = loc[0][qt][0] - loc[0][q][0];
            temp2 = loc[0][qt][1] - loc[0][q][1];
            temp = temp1 * temp1 + temp2 * temp2;
            temp = sqrt(temp);
                
            L[0][q].push_back(temp);
            L[1][q].push_back(0);

            Nx[0][q].push_back(- temp2 / temp);
            Ny[0][q].push_back(temp1 / temp);
            Nx[1][q].push_back(0);
            Ny[1][q].push_back(0);

            Pstar[q][0].push_back(0);
            Pstar[q][1].push_back(0);
        }
        
    }
    delete[] ntemp;
    delete[] ctemp;
    //划分量初始化完毕

    return;
}

void internalsolver(int t_l, int qm)
{
    int h, l, q, k, k2;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22, vk;
    
    //初始化矩阵
    double** A = new double* [2];
    double* B = new double [2];
    for (h=0; h<2; h++)
    {
        A[h] = new double [2];
    }
    for (h=0; h<2; h++)
    {
        for (l=0; l<2; l++)
        {
            A[h][l] = 0;
        }
        B[h] = 0;
    }


    //给矩阵赋值
    for (k=0; k< neighbornode[qm].size(); k++)
    {
        //先计算矩阵系数

        //k2 = k-1，注意是循环指标
        {
            k2 = k-1;
            if ( k == 0)
            {
                k2 = neighbornode[qm].size() - 1;
            }
        }

        //注意此时是内点所以每个邻域都有意义

        //计算alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            q = neighborcell[qm][k2];  //(i,j)点的第k个邻域的总体编号

            h = q / ny;
            l = q % ny;

            alpha2 = c[t_l][h][l] / tau[t_l][h][l];
            p2 = cellP[t_l][h][l];
            v21 = cellv[t_l][h][l][0];
            v22 = cellv[t_l][h][l][1];
        }
        
        //计算alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            q = neighborcell[qm][k];

            h = q / ny;
            l = q % ny;

            alpha1 = c[t_l][h][l] / tau[t_l][h][l];
            p1 = cellP[t_l][h][l];
            v11 = cellv[t_l][h][l][0];
            v12 = cellv[t_l][h][l][1];
        }

        //计算系数矩阵A
        {
            A[0][0] = A[0][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Nx[t_l][qm][k];
            A[0][1] = A[0][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][0] = A[1][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][1] = A[1][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Ny[t_l][qm][k] * Ny[t_l][qm][k];
        }

        //计算vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][qm][k] + v22 * Ny[t_l][qm][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][qm][k] + v12 * Ny[t_l][qm][k] );
            vk = vk / (alpha1 + alpha2);
        }

        //计算右端矩阵B
        {
            B[0] = B[0] + L[t_l][qm][k] * ( alpha1 + alpha2) * vk * Nx[t_l][qm][k];
            B[1] = B[1] + L[t_l][qm][k] * ( alpha1 + alpha2) * vk * Ny[t_l][qm][k];
        }

    }


    //解方程，二维矩阵直接求逆就行
    {
        double** Ainv = new double* [2];
        for (h=0; h<2; h++)
        {
            Ainv[h] = new double [2];
        }

        alpha1 = A[0][0] * A[1][1] - A[0][1] * A[1][0]; //alpha1为A的行列式

        Ainv[0][0] = A[1][1] / alpha1;
        Ainv[0][1] = - A[0][1] /alpha1;
        Ainv[1][0] = - A[1][0] / alpha1;
        Ainv[1][1] = A[0][0] / alpha1;

        for (h=0; h<2; h++)
        {
            nodev[qm][h] = 0;
            for (l=0; l<2; l++)
            {
                nodev[qm][h] = nodev[qm][h] + Ainv[h][l] * B[l];
            }
        }

        for (h=0; h<2; h++)
        {
            delete[] Ainv[h];
        }
        delete[] Ainv;
    }
    
    for (h=0; h<2; h++)
    {
        delete[] A[h];
    }
    delete[] A;
    delete[] B;
    return;
}

void getinternalPstar(int t_l, int qm)
{
    int i, j, h, l, q, k, k2;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22;

    if ( qm == 0)
    {
        i = 0;
        j = 0;
    }
    else{
        i = (q - 1) / (ny + 1) + 1;
        j = (q - 1) % (ny + 1);
    }

    for (k=0; k<neighbornode[qm].size(); k++)
    {
        //k2 = k-1，注意是循环指标
        {
            k2 = k-1;
            if ( k == 0)
            {
                k2 = neighbornode[qm].size() - 1;
            }
        }

        //计算alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            q = neighborcell[qm][k2];  //(i,j)点的第k个邻域的总体编号

            h = q / ny;
            l = q % ny;

            alpha2 = c[t_l][h][l] / tau[t_l][h][l];
            p2 = cellP[t_l][h][l];
            v21 = cellv[t_l][h][l][0];
            v22 = cellv[t_l][h][l][1];
        }
        
        //计算alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            q = neighborcell[qm][k];

            h = q / ny;
            l = q % ny;

            alpha1 = c[t_l][h][l] / tau[t_l][h][l];
            p1 = cellP[t_l][h][l];
            v11 = cellv[t_l][h][l][0];
            v12 = cellv[t_l][h][l][1];
        }

        //解P*_{q,k}^{k-1}
        {
            q = neighborcell[qm][k2];
            h = q / ny;
            l = q % ny;

            Pstar[qm][0][k] = 0;

            Pstar[qm][0][k] = Pstar[qm][0][k] + (nodev[qm][0] - cellv[t_l][h][l][0]) * Nx[t_l][qm][k];
            Pstar[qm][0][k] = Pstar[qm][0][k] + (nodev[qm][1] - cellv[t_l][h][l][1]) * Ny[t_l][qm][k];
            Pstar[qm][0][k] = - Pstar[qm][0][k] * alpha2 + p2;
        }

        //解P*_{q,k}^k
        {
            q = neighborcell[qm][k];
            h = q / ny;
            l = q % ny;

            Pstar[qm][1][k] = 0;

            Pstar[qm][1][k] = Pstar[qm][1][k] + (nodev[qm][0] - cellv[t_l][h][l][0]) * Nx[t_l][qm][k];
            Pstar[qm][1][k] = Pstar[qm][1][k] + (nodev[qm][1] - cellv[t_l][h][l][1]) * Ny[t_l][qm][k];
            Pstar[qm][1][k] = Pstar[qm][1][k] * alpha1 + p1;
        }
    }

    return;
}

void boundarysolverNP(int t_l, int qm)    //此时算法为给定边界外法向速度的解法
{
    int h, l, k;
    double alpha1;
    
    //初始化矩阵
    double** A = new double* [2];
    double* B = new double [2];
    for (h=0; h<2; h++)
    {
        A[h] = new double [2];
    }
    for (h=0; h<2; h++)
    {
        for (l=0; l<2; l++)
        {
            A[h][l] = 0;
        }
        B[h] = 0;
    }

    //给矩阵赋值
    //计算系数矩阵A
    {
        k = neighbornode[qm].size() - 1;
        A[0][0] = - Nx[t_l][qm][0];
        A[0][1] = - Ny[t_l][qm][0];
        A[1][0] = Nx[t_l][qm][k]; //  和论文算法差一个-号
        A[1][1] = Ny[t_l][qm][k]; //和论文算法差一个-号
    }

    //计算右端矩阵B
    {
        B[0] = vn;
        B[1] = vn;
    }

    //解方程，二维矩阵直接求逆就行
    {
        double** Ainv = new double* [2];
        for (h=0; h<2; h++)
        {
            Ainv[h] = new double [2];
        }

        alpha1 = A[0][0] * A[1][1] - A[0][1] * A[1][0]; //alpha1为A的行列式

        Ainv[0][0] = A[1][1] / alpha1;
        Ainv[0][1] = - A[0][1] /alpha1;
        Ainv[1][0] = - A[1][0] / alpha1;
        Ainv[1][1] = A[0][0] / alpha1;

        for (h=0; h<2; h++)
        {
            nodev[qm][h] = 0;
            for (l=0; l<2; l++)
            {
                nodev[qm][h] = nodev[qm][h] + Ainv[h][l] * B[l];
            }
        }

        for (h=0; h<2; h++)
        {
            delete[] Ainv[h];
        }
        delete[] Ainv;
    }
    

    for (h=0; h<2; h++)
    {
        delete[] A[h];
    }
    delete[] A;
    delete[] B;

    return;
}

void boundarysolverP(int t_l, int qm)
{
    int h, l, q, k, k2;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22, vk;
    
    //初始化矩阵
    double** A = new double* [3];
    double* B = new double [3];
    for (h=0; h<3; h++)
    {
        A[h] = new double [3];
    }
    for (h=0; h<3; h++)
    {
        for (l=0; l<3; l++)
        {
            A[h][l] = 0;
        }
        B[h] = 0;
    }

    //给矩阵赋值
    for (k=0; k< neighbornode[qm].size(); k++)
    {
        //先计算矩阵系数

        //k2 = k-1，注意此处不是循环指标
        k2 = k - 1;

        //计算alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            if (  k == 0 )  //-1邻点对应的网格为虚网格
            {
                alpha2 = 0;
                p2 = 0;
                v21 = 0;
                v22 = 0;
            }
            else{
                q = neighborcell[qm][k2];  //(i,j)点的第k个邻域的总体编号

                h = q / ny;
                l = q % ny;

                alpha2 = c[t_l][h][l] / tau[t_l][h][l];
                p2 = cellP[t_l][h][l];
                v21 = cellv[t_l][h][l][0];
                v22 = cellv[t_l][h][l][1];
            }
            
        }
        
        //计算alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            if ( k == neighbornode[qm].size() - 1 ) //第k个邻点对应的网格为虚网格
            {
                alpha1 = 0;
                p1 = 0;
                v11 = 0;
                v12 = 0;
            }
            else{
                q = neighborcell[qm][k];

                h = q / ny;
                l = q % ny;

                alpha1 = c[t_l][h][l] / tau[t_l][h][l];
                p1 = cellP[t_l][h][l];
                v11 = cellv[t_l][h][l][0];
                v12 = cellv[t_l][h][l][1];
            }
            
        }

        //计算系数矩阵A
        {
            A[0][0] = A[0][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Nx[t_l][qm][k];
            A[0][1] = A[0][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][0] = A[1][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][1] = A[1][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Ny[t_l][qm][k] * Ny[t_l][qm][k];
            //先只给和k有关的赋值
        }

        //计算vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][qm][k] + v22 * Ny[t_l][qm][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][qm][k] + v12 * Ny[t_l][qm][k] );
        }

        //计算右端矩阵B
        {
            B[0] = B[0] + L[t_l][qm][k] * vk * Nx[t_l][qm][k];
            B[1] = B[1] + L[t_l][qm][k] * vk * Ny[t_l][qm][k];
        }

    }
    
    //给剩下的赋值
    {
        k2 = neighbornode[qm].size() - 1;

        A[0][2] = - L[t_l][qm][0] * Nx[t_l][qm][0] + L[t_l][qm][k2] * Nx[t_l][qm][k2];
        A[2][0] = A[0][2];
        A[1][2] = - L[t_l][qm][0] * Ny[t_l][qm][0] + L[t_l][qm][k2] * Ny[t_l][qm][k2];
        A[2][1] = A[1][2];
        A[2][2] = 0;

        B[2] = L[t_l][qm][0] * vn + L[t_l][qm][k2] * vn;
    }

    //最后解线性方程组，注意边界外法向量共线时解是三维的，只取速度的部分

    {
        for (h =1; h<3; h++)
        {
            A[h][0]=A[h][0]/A[0][0];
        }
        for (h=1; h<3; h++)
        {
            for (l=h; l<3; l++)
            {
                for (k=0; k<h; k++)
                {
                    A[h][l] = A[h][l] - A[h][k] * A[k][l];
                }
            }
            for (l=h+1; l<3; l++)
            {
                for (k=0; k<h; k++)
                {
                    A[l][h] = A[l][h] - A[l][k] * A[k][h];
                }
                    A[l][h] = A[l][h] / A[h][h];
            }

            for (l = 0; l< h; l++)
            {
                B[h] = B[h] - A[h][l] * B[l];
            }
        }

        B[2] = B[2] / A[2][2];
        for (h = 1; h>=0; h--)
        {
            for (l=2; l>h; l--)
            {
                B[h] = B[h] - A[h][l] * B[l];
            }
            B[h] = B[h] / A[h][h];
        }
    }
    
    nodev[qm][0] = B[0];
    nodev[qm][1] = B[1];

    for (h=0; h<3; h++)
    {
        delete[] A[h];
    }
    delete[] A;
    delete[] B;
}

void anglesolver1(int t_l, int qm)
{
    int h, l, k, q;
    double alpha, temp1, temp2;
    double** A = new double* [3];
    double* B = new double [3];
    for (h=0; h<3; h++)
    {
        A[h] = new double [3];
    }
    for (h=0; h<3; h++)
    {
        for (l=0; l<3; l++)
        {
            A[h][l] = 0;
        }
        B[h] = 0;
    }

    
    q = neighborcell[qm][0];
    h = q / ny;
    l = q % ny;
    alpha = c[t_l][h][l] / tau[t_l][h][l];

    temp1 = cellv[t_l][h][l][0] * Nx[t_l][qm][0] + cellv[t_l][h][l][1] * Ny[t_l][qm][0];
    temp2 = cellv[t_l][h][l][0] * Nx[t_l][qm][1] + cellv[t_l][h][l][1] * Ny[t_l][qm][1];

    
    //给A赋值
    {
        A[0][0] = L[t_l][qm][0] * alpha * (Nx[t_l][qm][0] * Nx[t_l][qm][0]);
        A[0][0] = A[0][0] + L[t_l][qm][1] * alpha * (Nx[t_l][qm][1] * Nx[t_l][qm][1]);

        A[0][1] = L[t_l][qm][0] * alpha * (Nx[t_l][qm][0] * Ny[t_l][qm][0]);
        A[0][1] = A[0][1] + L[t_l][qm][1] * alpha * (Nx[t_l][qm][1] * Ny[t_l][qm][1]);

        A[0][2] = L[t_l][qm][1] * Nx[t_l][qm][1];

        A[1][0] = A[0][1];

        A[1][1] = L[t_l][qm][0] * alpha * (Ny[t_l][qm][0] * Ny[t_l][qm][0]);
        A[1][1] = A[1][1] + L[t_l][qm][1] * alpha * (Ny[t_l][qm][1] * Ny[t_l][qm][1]);

        A[1][2] = L[t_l][qm][1] * Ny[t_l][qm][1];

        A[2][0] = Nx[t_l][qm][1];
        A[2][1] = Ny[t_l][qm][1];
        A[2][2] = 0;
    }

    //给B赋值
    {
        B[0] = - L[t_l][qm][0] * cellP[t_l][h][l] * Nx[t_l][qm][0];
        B[0] = B[0] + L[t_l][qm][0] * alpha * temp1 * Nx[t_l][qm][0];
        B[0] = B[0] + L[t_l][qm][1] * cellP[t_l][h][l] * Nx[t_l][qm][1];
        B[0] = B[0] + L[t_l][qm][1] * alpha * temp2 * Nx[t_l][qm][1];
        B[0] = B[0] + L[t_l][qm][0] * P_0 * Nx[t_l][qm][0];

        B[1] = - L[t_l][qm][0] * cellP[t_l][h][l] * Ny[t_l][qm][0];
        B[1] = B[1] + L[t_l][qm][0] * alpha * temp1 * Ny[t_l][qm][0];
        B[1] = B[1] + L[t_l][qm][1] * cellP[t_l][h][l] * Ny[t_l][qm][1];
        B[1] = B[1] + L[t_l][qm][1] * alpha * temp2 * Ny[t_l][qm][1];
        B[1] = B[1] + L[t_l][qm][0] * P_0 * Ny[t_l][qm][0];

        B[2] = vn;
    }

    //解方程
    {
        for (h =1; h<3; h++)
        {
            A[h][0]=A[h][0]/A[0][0];
        }
        for (h=1; h<3; h++)
        {
            for (l=h; l<3; l++)
            {
                for (k=0; k<h; k++)
                {
                    A[h][l] = A[h][l] - A[h][k] * A[k][l];
                }
            }
            for (l=h+1; l<3; l++)
            {
                for (k=0; k<h; k++)
                {
                    A[l][h] = A[l][h] - A[l][k] * A[k][h];
                }
                    A[l][h] = A[l][h] / A[h][h];
            }

            for (l = 0; l< h; l++)
            {
                B[h] = B[h] - A[h][l] * B[l];
            }
        }

        B[2] = B[2] / A[2][2];
        for (h = 1; h>=0; h--)
        {
            for (l=2; l>h; l--)
            {
                B[h] = B[h] - A[h][l] * B[l];
            }
            B[h] = B[h] / A[h][h];
        }
    }
    
    nodev[qm][0] = B[0];
    nodev[qm][1] = B[1];

    for (h=0; h<3; h++)
    {
        delete[] A[h];
    }
    delete[] A;
    delete[] B;
    return;
}

void anglesolver2(int t_l, int qm)
{
    int h, l, k, q;
    double alpha, temp1, temp2;
    double** A = new double* [3];
    double* B = new double [3];
    for (h=0; h<3; h++)
    {
        A[h] = new double [3];
    }
    for (h=0; h<3; h++)
    {
        for (l=0; l<3; l++)
        {
            A[h][l] = 0;
        }
        B[h] = 0;
    }

    
    q = neighborcell[qm][0];
    h = q / ny;
    l = q % ny;
    alpha = c[t_l][h][l] / tau[t_l][h][l];

    temp1 = cellv[t_l][h][l][0] * Nx[t_l][qm][0] + cellv[t_l][h][l][1] * Ny[t_l][qm][0];
    temp2 = cellv[t_l][h][l][0] * Nx[t_l][qm][1] + cellv[t_l][h][l][1] * Ny[t_l][qm][1];

    
    //给A赋值
    {
        A[0][0] = L[t_l][qm][0] * alpha * (Nx[t_l][qm][0] * Nx[t_l][qm][0]);
        A[0][0] = A[0][0] + L[t_l][qm][1] * alpha * (Nx[t_l][qm][1] * Nx[t_l][qm][1]);

        A[0][1] = L[t_l][qm][0] * alpha * (Nx[t_l][qm][0] * Ny[t_l][qm][0]);
        A[0][1] = A[0][1] + L[t_l][qm][1] * alpha * (Nx[t_l][qm][1] * Ny[t_l][qm][1]);

        A[0][2] = - L[t_l][qm][0] * Nx[t_l][qm][0];

        A[1][0] = A[0][1];

        A[1][1] = L[t_l][qm][0] * alpha * (Ny[t_l][qm][0] * Ny[t_l][qm][0]);
        A[1][1] = A[1][1] + L[t_l][qm][1] * alpha * (Ny[t_l][qm][1] * Ny[t_l][qm][1]);

        A[1][2] = - L[t_l][qm][0] * Ny[t_l][qm][0];

        A[2][0] = - Nx[t_l][qm][0];
        A[2][1] = - Ny[t_l][qm][0];
        A[2][2] = 0;
    }

    //给B赋值
    {
        B[0] = - L[t_l][qm][0] * cellP[t_l][h][l] * Nx[t_l][qm][0];
        B[0] = B[0] + L[t_l][qm][0] * alpha * temp1 * Nx[t_l][qm][0];
        B[0] = B[0] + L[t_l][qm][1] * cellP[t_l][h][l] * Nx[t_l][qm][1];
        B[0] = B[0] + L[t_l][qm][1] * alpha * temp2 * Nx[t_l][qm][1];
        B[0] = B[0] - L[t_l][qm][1] * P_0 * Nx[t_l][qm][1];

        B[1] = - L[t_l][qm][0] * cellP[t_l][h][l] * Ny[t_l][qm][0];
        B[1] = B[1] + L[t_l][qm][0] * alpha * temp1 * Ny[t_l][qm][0];
        B[1] = B[1] + L[t_l][qm][1] * cellP[t_l][h][l] * Ny[t_l][qm][1];
        B[1] = B[1] + L[t_l][qm][1] * alpha * temp2 * Ny[t_l][qm][1];
        B[1] = B[1] - L[t_l][qm][1] * P_0 * Ny[t_l][qm][1];

        B[2] = vn;
    }

    //解方程
    {
        for (h =1; h<3; h++)
        {
            A[h][0]=A[h][0]/A[0][0];
        }
        for (h=1; h<3; h++)
        {
            for (l=h; l<3; l++)
            {
                for (k=0; k<h; k++)
                {
                    A[h][l] = A[h][l] - A[h][k] * A[k][l];
                }
            }
            for (l=h+1; l<3; l++)
            {
                for (k=0; k<h; k++)
                {
                    A[l][h] = A[l][h] - A[l][k] * A[k][h];
                }
                    A[l][h] = A[l][h] / A[h][h];
            }

            for (l = 0; l< h; l++)
            {
                B[h] = B[h] - A[h][l] * B[l];
            }
        }

        B[2] = B[2] / A[2][2];
        for (h = 1; h>=0; h--)
        {
            for (l=2; l>h; l--)
            {
                B[h] = B[h] - A[h][l] * B[l];
            }
            B[h] = B[h] / A[h][h];
        }
    }
    
    nodev[qm][0] = B[0];
    nodev[qm][1] = B[1];

    for (h=0; h<3; h++)
    {
        delete[] A[h];
    }
    delete[] A;
    delete[] B;
    return;
}

void boundarywithP(int t_l, int qm)
{
    int h, l, k, k2, q;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22, vk;
    double** A = new double* [2];
    double* B = new double [2];
    for (h=0; h<2; h++)
    {
        A[h] = new double [2];
    }
    for (h=0; h<2; h++)
    {
        for (l=0; l<2; l++)
        {
            A[h][l] = 0;
        }
        B[h] = 0;
    }

    //先给矩阵赋值
    for (k=0; k< neighbornode[qm].size(); k++)
    {
        //先计算矩阵系数

        //k2 = k-1，注意此处不是循环指标
        k2 = k - 1;

        //计算alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            if (  k == 0 )  //-1邻点对应的网格为虚网格
            {
                alpha2 = 0;
                p2 = P_0;
                v21 = 0;
                v22 = 0;
            }
            else{
                q = neighborcell[qm][k2];  //(i,j)点的第k个邻域的总体编号

                h = q / ny;
                l = q % ny;

                alpha2 = c[t_l][h][l] / tau[t_l][h][l];
                p2 = cellP[t_l][h][l];
                v21 = cellv[t_l][h][l][0];
                v22 = cellv[t_l][h][l][1];
            }
            
        }
        
        //计算alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            if ( k == neighbornode[qm].size() - 1 ) //第k个邻点对应的网格为虚网格
            {
                alpha1 = 0;
                p1 = P_0;
                v11 = 0;
                v12 = 0;
            }
            else{
                q = neighborcell[qm][k];

                h = q / ny;
                l = q % ny;

                alpha1 = c[t_l][h][l] / tau[t_l][h][l];
                p1 = cellP[t_l][h][l];
                v11 = cellv[t_l][h][l][0];
                v12 = cellv[t_l][h][l][1];
            }
            
        }

        //计算系数矩阵A
        {
            A[0][0] = A[0][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Nx[t_l][qm][k];
            A[0][1] = A[0][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][0] = A[1][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][1] = A[1][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Ny[t_l][qm][k] * Ny[t_l][qm][k];
            //先只给和k有关的赋值
        }

        //计算vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][qm][k] + v22 * Ny[t_l][qm][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][qm][k] + v12 * Ny[t_l][qm][k] );
        }

        //计算右端矩阵B
        {
            B[0] = B[0] + L[t_l][qm][k] * vk * Nx[t_l][qm][k];
            B[1] = B[1] + L[t_l][qm][k] * vk * Ny[t_l][qm][k];
        }

    }

     //解方程，二维矩阵直接求逆就行
    {
        double** Ainv = new double* [2];
        for (h=0; h<2; h++)
        {
            Ainv[h] = new double [2];
        }

        alpha1 = A[0][0] * A[1][1] - A[0][1] * A[1][0]; //alpha1为A的行列式

        Ainv[0][0] = A[1][1] / alpha1;
        Ainv[0][1] = - A[0][1] /alpha1;
        Ainv[1][0] = - A[1][0] / alpha1;
        Ainv[1][1] = A[0][0] / alpha1;

        for (h=0; h<2; h++)
        {
            nodev[qm][h] = 0;
            for (l=0; l<2; l++)
            {
                nodev[qm][h] = nodev[qm][h] + Ainv[h][l] * B[l];
            }
        }

        for (h=0; h<2; h++)
        {
            delete[] Ainv[h];
        }
        delete[] Ainv;
    }

    for (h=0; h<2; h++)
    {
        delete[] A[h];
    }
    delete[] A;
    delete[] B;
    return;
}

void boundarysolver2(int t_l, int qm)
{
    
    if ( qm == (nx - 1)*(ny+1) + 1 )    //q为(nx,0)的情况
    {
        anglesolver1(t_l, qm);
    }
    else if ( qm == n ) //(nx, ny)的情况
    {
        anglesolver2(t_l, qm);
    }
    else{
        boundarywithP(t_l, qm);
    }

    return;
}

void getboundaryPstar(int t_l, int qm)
{
    int h, l, q, k, k2;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22;

    for (k=0; k<neighbornode[qm].size(); k++)
    {
        //k2 = k-1
        k2 = k - 1;

        //计算alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            if ( k == 0 )
            {
                alpha2 = 0;
                p2 = 0;
                v21 = 0;
                v22 = 0;
            }
            else{
                q = neighborcell[qm][k2];  //(i,j)点的第k个邻域的总体编号

                h = q / ny;
                l = q % ny;

                alpha2 = c[t_l][h][l] / tau[t_l][h][l];
                p2 = cellP[t_l][h][l];
                v21 = cellv[t_l][h][l][0];
                v22 = cellv[t_l][h][l][1];
            }
        }
        
        //计算alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            if ( k == neighbornode[qm].size() - 1 )
            {
                alpha1 = 0;
                p1 = 0;
                v11 = 0;
                v12 = 0;
            }
            else{
                q = neighborcell[qm][k];

                h = q / ny;
                l = q % ny;

                alpha1 = c[t_l][h][l] / tau[t_l][h][l];
                p1 = cellP[t_l][h][l];
                v11 = cellv[t_l][h][l][0];
                v12 = cellv[t_l][h][l][1];
            }
        }

        //解P*_{q,k}^{k-1}
        {
            if ( k == 0 )
            {
                Pstar[qm][0][k] = 0;
            }
            else{
                q = neighborcell[qm][k2];
                h = q / ny;
                l = q % ny;

                Pstar[qm][0][k] = 0;

                Pstar[qm][0][k] = Pstar[qm][0][k] + (nodev[qm][0] - cellv[t_l][h][l][0]) * Nx[t_l][qm][k];
                Pstar[qm][0][k] = Pstar[qm][0][k] + (nodev[qm][1] - cellv[t_l][h][l][1]) * Ny[t_l][qm][k];
                Pstar[qm][0][k] = - Pstar[qm][0][k] * alpha2 + p2;
            }
        }

        //解P*_{q,k}^k
        {
            if ( k == neighbornode[qm].size() - 1 )
            {
                Pstar[qm][1][k] = 0;
            }
            else{
                q = neighborcell[qm][k];
                h = q / ny;
                l = q % ny;

                Pstar[qm][1][k] = 0;

                Pstar[qm][1][k] = Pstar[qm][1][k] + (nodev[qm][0] - cellv[t_l][h][l][0]) * Nx[t_l][qm][k];
                Pstar[qm][1][k] = Pstar[qm][1][k] + (nodev[qm][1] - cellv[t_l][h][l][1]) * Ny[t_l][qm][k];
                Pstar[qm][1][k] = Pstar[qm][1][k] * alpha1 + p1;
            }
        }
    }

    return;
}

void nodalsolver(int t_l)
{
    int i, j, h, l, q, qm;
    double x1, x2, y1, y2, temp;

    for (qm=0; qm<=n; qm++)
    {
        if ( qm == 0)
        {
            i = 0;
            j = 0;
        }
        else{
            i = (qm - 1) / (ny + 1) + 1;
            j = (qm - 1) % (ny + 1);
        }

        //判断是否是边界点，根据不同点计算
        if (boundary[qm] == 0)    //是内点
        {
            internalsolver(t_l,qm);    //计算第(i,j)个点的速度
            getinternalPstar(t_l,qm);  //计算第(i,j)个点的边上的压力
        }
        else if(boundary[qm] == 1)  //边界条件由法向速度给出
        {
            //判断边界是否平行
            double e = 1e-10;

            q = neighbornode[qm].size()-1;
                
            x1 = Nx[t_l][qm][0];
            y1 = Ny[t_l][qm][0];
            x2 = Nx[t_l][qm][q];
            y2 = Ny[t_l][qm][q];

            temp = x1 * y2 - x2 * y1;
            temp = abs(temp);

            if (temp >e)    //不平行时
            {
                boundarysolverNP(t_l,qm);
            }
            else{
                boundarysolverP(t_l,qm);
            }
            getboundaryPstar(t_l,qm);
        }
        else{   //边界条件同时由法向速度和压力给出
            boundarysolver2(t_l,qm);
            getboundaryPstar(t_l,qm);
        }

    }
    return;
}

double choosedt(int t_l, double dtn)
{
    double dt, te, tv, lambda, temp;
    int h, l, k, k2, qm, qm2;
    
    //CFL criterion
    {
        te = 1000;

        int* vertex = new int [4];

        for (h=0; h<nx; h++)
        {
            for (l=0; l<ny; l++)
            {
                //先记录(h,l)网格每个顶点的编号
                {
                    if (h == 0)
                    {
                        vertex[0] = 0;
                        vertex[1] = l + 1;
                        vertex[2] = l + 2;
                        vertex[3] = 0;
                    }
                    else{
                        vertex[0] = (h - 1) * (ny + 1) + l + 1;
                        vertex[1] = (h + 1 - 1) * (ny + 1) + l + 1;
                        vertex[2] = (h + 1 - 1) * (ny + 1) + l + 1 + 1;
                        vertex[3] = (h - 1) * (ny + 1) + l + 1 + 1;
                    }
                    
                }

                //选出每个网格中最小的边长
                {
                    lambda = 10000;
                    for (k=0; k<4; k++)
                    {
                        //qm为(h,l)的第k个顶点
                        qm = vertex[k];
                        
                        //计算其他顶点到(i,j)点的最短距离
                        k2 = k + 1;
                        if ( k2 > 3)
                        {
                            k2 = k2 - 4;
                        }

                        do{
                            qm2 = vertex[k2];

                            temp = (loc[t_l][qm][0] - loc[t_l][qm2][0]) * (loc[t_l][qm][0] - loc[t_l][qm2][0]);
                            temp = temp + (loc[t_l][qm][1] - loc[t_l][qm2][1]) * (loc[t_l][qm][1] - loc[t_l][qm2][1]);
                            temp = sqrt(temp);

                            if (temp < lambda && temp != 0)
                            {
                                lambda = temp;
                            }
                            k2++;
                            if ( k2 > 3)
                            {
                                k2 = k2 - 4;
                            }

                        }while(k2 != k);
                    }
                }

                //选取te
                lambda = lambda / c[t_l][h][l];

                if (lambda < te)
                {
                    te = lambda;
                }
            }
        }
        
        te = Ce * te;
        
        delete[] vertex;

    }

    //criterion on the variation of volume
    {
        tv = 1000;

        int* vertex = new int [4];

        for (h=0; h<nx; h++)
        {
            for (l=0; l<ny; l++)
            {
                //先记录(h,l)网格每个顶点的编号
                {
                    if (h == 0)
                    {
                        vertex[0] = 0;
                        vertex[1] = l + 1;
                        vertex[2] = l + 2;
                        vertex[3] = 0;
                    }
                    else{
                        vertex[0] = (h - 1) * (ny + 1) + l + 1;
                        vertex[1] = (h + 1 - 1) * (ny + 1) + l + 1;
                        vertex[2] = (h + 1 - 1) * (ny + 1) + l + 1 + 1;
                        vertex[3] = (h - 1) * (ny + 1) + l + 1 + 1;
                    }
                }

                //计算每个区域的面积变化率
                {
                    temp = 0;   //temp记录面积变化率
                    for (k=0; k<4; k++)
                    {
                        //qm为(h,l)的第k个顶点
                        qm = vertex[k];

                        //k2 = k + 1, 注意是循环指标
                        k2 = k + 1;
                        if ( k2 > 3)
                        {
                            k2 = k2 - 4;
                        }

                        qm2 = vertex[k2];

                        lambda = nodev[qm][0] * loc[t_l][qm2][1];
                        lambda = lambda + nodev[qm2][1] * loc[t_l][qm][0];
                        lambda = lambda - nodev[qm2][0] * loc[t_l][qm][1];
                        lambda = lambda - nodev[qm][1] * loc[t_l][qm2][0];

                        temp = temp + lambda;
                        
                    }
                    temp = temp * 0.5;
                    temp = abs(temp);
                }

                //选取tv
                
                temp = A[t_l][h][l] / temp;

                if (temp < tv && temp != 0 )
                {
                    tv = temp;
                }
            }
        }

        tv = Cv * tv;
        
        delete[] vertex;

    }

    //选取dt
    {
        dt = Cm * dtn;
        
        if (te < dt)
        {
            dt = te;
        }

        if (tv < dt)
        {
            dt = tv;
        }
    }
    
    return dt;
}

void updategeometric(int t_l, int t_n, double dt)
{
    int h, l, k, k2, qm, qm2;
    double temp;

    //先更新节点位置
    {
        for (qm=0; qm<= n; qm++)
        {
            loc[t_n][qm][0] = loc[t_l][qm][0] + dt * nodev[qm][0];
            loc[t_n][qm][1] = loc[t_l][qm][1] + dt * nodev[qm][1];
        }
    }

    //再更新网格面积
    {
        int* vertex = new int [4];

        for (h=0; h<nx; h++)
        {
            for (l=0; l <ny; l++)
            {
                A[t_n][h][l] = 0;
                //先记录(h,l)网格每个顶点的编号
                {
                    if (h == 0)
                    {
                        vertex[0] = 0;
                        vertex[1] = l + 1;
                        vertex[2] = l + 2;
                        vertex[3] = 0;
                    }
                    else{
                        vertex[0] = (h - 1) * (ny + 1) + l + 1;
                        vertex[1] = (h + 1 - 1) * (ny + 1) + l + 1;
                        vertex[2] = (h + 1 - 1) * (ny + 1) + l + 1 + 1;
                        vertex[3] = (h - 1) * (ny + 1) + l + 1 + 1;
                    }
                }
                
                for (k=0; k < 3; k++)
                {
                    //k2 = k+1,注意是循环指标
                    k2 = k + 1;
                    if (k2 > 3)
                    {
                        k2 = k2 - 4;
                    }

                    qm = vertex[k];

                    qm2 = vertex[k2];

                    if (qm == 0 && qm2 == 0)
                    {
                        temp = 0;
                    }
                    else{
                        temp = ( loc[t_l][qm][0] + dt * nodev[qm][0] ) * ( loc[t_l][qm2][1] + dt * nodev[qm2][1] );
                        temp = temp - ( loc[t_l][qm][1] + dt * nodev[qm][1] ) * ( loc[t_l][qm2][0] + dt * nodev[qm2][0] );
                    }
                   
                    A[t_n][h][l] = A[t_l][h][l] + temp * 0.5;

                    if (A[t_n][h][l] < 0)
                    {
                        cout<<endl<<"\t"<<"\t"<<h<<"\t"<<"\t"<<l<<endl<<endl;
                    }
                }
            }
        }

        delete[] vertex;
    }

    //再更新区域的边长和外法向量
    {
        
        for (qm=0; qm <= n; qm++)
        {
            for (k=0; k<neighbornode[qm].size(); k++)
            {
                qm2 = neighbornode[qm][k];

                //先更新边长
                temp = ( loc[t_n][qm2][0] - loc[t_n][qm][0] ) * ( loc[t_n][qm2][0] - loc[t_n][qm][0] );
                temp = temp + ( loc[t_n][qm2][1] - loc[t_n][qm][1] ) * ( loc[t_n][qm2][1] - loc[t_n][qm][1] );
                temp = sqrt(temp);

                L[t_n][qm][k] = temp;

                //再更新外法向量
                Nx[t_n][qm][k] = - (loc[t_n][qm2][1] - loc[t_n][qm][1]) / temp;
                Ny[t_n][qm][k] = (loc[t_n][qm2][0] - loc[t_n][qm][0]) / temp;
            }
        }
    }
    return ;
}

void updatephysicalvariables(int t_l, int t_n, double dt)
{
    int h, l, k, k2, r, qm, qm2;
    double temp, temp1, temp2;

    int* vertex = new int [4];

    for (h=0; h < nx; h++)
    {
        for (l=0; l<ny; l++)
        {
            //记录顶点序号
            {
                if (h == 0)
                {
                    vertex[0] = 0;
                    vertex[1] = l + 1;
                    vertex[2] = l + 2;
                    vertex[3] = 0;
                }
                else{
                    vertex[0] = (h - 1) * (ny + 1) + l + 1;
                    vertex[1] = (h + 1 - 1) * (ny + 1) + l + 1;
                    vertex[2] = (h + 1 - 1) * (ny + 1) + l + 1 + 1;
                    vertex[3] = (h - 1) * (ny + 1) + l + 1 + 1;
                }
            }

            //先更新比容
            {
                temp = 0;
                for (k=0; k< 4; k++)
                {
                    qm = vertex[k];   //qm为(h,l)的第k个顶点

                    k2 = k + 1;
                    if (k2 > 3)
                    {
                        k2 = k2 - 4;
                    }

                    qm2 = vertex[k2];  //qm2为(h,l)的第k+1个顶点

                    //确认第k条边长和法向量
                    for (r=0; r<neighbornode[qm].size(); r++)
                    {
                        if ( neighbornode[qm][r] == vertex[k2] )  //确认k+1顶点是k顶点的第几个邻点
                        {
                            temp1 = Nx[t_l][qm][r] * (nodev[qm][0] + nodev[qm2][0]) + Ny[t_l][qm][r] * (nodev[qm][1] + nodev[qm2][1]);
                            temp1 = temp1 * L[t_l][qm][r];

                            temp2 = Nx[t_n][qm][r] * (nodev[qm][0] + nodev[qm2][0]) + Ny[t_n][qm][r] * (nodev[qm][1] + nodev[qm2][1]);
                            temp2 = temp2 * L[t_n][qm][r];

                            temp = temp + temp1 + temp2;

                            break;
                        }
                        else{
                            ;
                        }
                    }
                    
                }
                
                temp = - temp * dt * 0.25;  //负号由区域外法向量和N的方向相反得到
                temp = temp / m[h][l];
                
                temp = temp + tau[t_l][h][l];

                tau[t_n][h][l] = temp;
            }

            //再更新速度
            {
                temp1 = 0;
                temp2 = 0;

                for (k=0; k< 4; k++)
                {
                    qm = vertex[k];

                    //先算L_{r,r+1}^n * P_{r,r+1/2}^{*,i} * N_{r,r+1}^n
                    {
                        k2 = k + 1;
                        if (k2 > 3)
                        {
                            k2 = k2 - 4;
                        }

                        for (r=0; r < neighbornode[qm].size(); r++)
                        {
                            if ( neighbornode[qm][r] == vertex[k2] )
                            {
                                //注意区域外法向量和N的方向关系
                                temp1 = temp1 - L[t_l][qm][r] * Pstar[qm][1][r] * Nx[t_l][qm][r];
                                temp2 = temp2 - L[t_l][qm][r] * Pstar[qm][1][r] * Ny[t_l][qm][r];
                                break;
                            }
                        }
                    }
                    
                    //再算L_{r,r+1}^n * P_{r+1/2, r+1}^{*,i} * N_{r,r+1}^n
                    {
                        k2 = k - 1;
                        if (k2 < 0)
                        {
                            k2 = k2 + 4;
                        }

                        for (r=0; r<neighbornode[qm].size(); r++)
                        {
                            if ( neighbornode[qm][r] == vertex[k2] )
                            {
                                temp1 = temp1 + L[t_l][qm][r] * Pstar[qm][0][r] * Nx[t_l][qm][r];
                                temp2 = temp2 + L[t_l][qm][r] * Pstar[qm][0][r] * Ny[t_l][qm][r];
                                break;
                            }
                        }
                    }
                    
                }

                temp1 = temp1 * dt * 0.5;
                temp2 = temp2 * dt * 0.5;

                temp1 = - temp1 / m[h][l];
                temp2 = - temp2 / m[h][l];

                temp1 = temp1 + cellv[t_l][h][l][0];
                temp2 = temp2 + cellv[t_l][h][l][1];

                cellv[t_n][h][l][0] = temp1;
                cellv[t_n][h][l][1] = temp2;
                
            }
            
            //再更新能量
            {
                temp = 0;

                for (k=0; k<4; k++)
                {
                    qm = vertex[k];

                    //先计算L_{r,r+1}^n * P_{r,r+1/2}^{*,i} * V*_r * N_{r,r+1}^n
                    {
                        k2 = k + 1;
                        if (k2 > 3)
                        {
                            k2 = k2 - 4;
                        }

                        for (r=0; r<neighbornode[qm].size(); r++)
                        {
                            if ( neighbornode[qm][r] == vertex[k2] )
                            {
                                temp1 = nodev[qm][0] * Nx[t_l][qm][r] + nodev[qm][1] * Ny[t_l][qm][r];
                                temp = temp - L[t_l][qm][r] * Pstar[qm][1][r] * temp1;  //注意是区域的外法向量和N的方向关系
                                break;
                            }
                        }

                    }

                    //再算L_{r,r+1}^n * P_{r+1/2, r+1}^{*,i} * V*_{r+1} * N_{r,r+1}^n
                    //注意循环指标的变化
                    {
                        k2 = k - 1;
                        if (k2 < 0)
                        {
                            k2 = k2 + 4;
                        }

                        for (r=0; r<neighbornode[qm].size(); r++)
                        {
                            if ( neighbornode[qm][r] == vertex[k2] )
                            {
                                temp1 = nodev[qm][0] * Nx[t_l][qm][r] + nodev[qm][1] * Ny[t_l][qm][r];
                                temp = temp + L[t_l][qm][r] * Pstar[qm][0][r] * temp1;
                                break;
                            }
                        }
                    }


                }
                temp = - temp * dt * 0.5;
                temp = temp / m[h][l];
                temp = temp + E[t_l][h][l];
                E[t_n][h][l] = temp;
            }
        }
    }

    delete[] vertex;

    return;
}

void getstates(int t_n)
{
    int h, l;
    double temp;

    for (h=0; h<nx; h++)
    {
        for (l=0; l< ny; l++)
        {
            //先算每个单元的内能
            {
                temp = cellv[t_n][h][l][0] * cellv[t_n][h][l][0] + cellv[t_n][h][l][1] * cellv[t_n][h][l][1];
                e[t_n][h][l] = E[t_n][h][l] - 0.5 * temp;
            }

            //再根据状态方程求单元压力
            {
                cellP[t_n][h][l] = ( gamma[h][l] - 1 ) * e[t_n][h][l];
                cellP[t_n][h][l] = cellP[t_n][h][l] / tau[t_n][h][l]; 
            }

            //求解声速
            {
                c[t_n][h][l] = (gamma[h][l] - 1) * e[t_n][h][l];
                c[t_n][h][l] = sqrt(c[t_n][h][l]);
            }
        }
    }

    return ;
}
