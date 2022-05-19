#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "initial.h"

using namespace std;

//************���涨����Ҫ�õ�������***********//

//���ڽڵ����

static double loc[2][n+1][2]; //loc[a][qm][k]Ϊ��(i,j)���ڵ�taʱ�̵ĵ�k��������ta=tn��tn+1
static double nodev[n+1][2];    //nodev[qm][k]Ϊ��(i,j)����ĵ�k���ٶȷ����������㷨ֻ��Ҫһ������洢
static int boundary[n+1];   //����0�������ڵ㣬1�Ǳ��ϱ߽�㣬2�ǻ��ϱ߽��

//�����������

static double tau[2][nx][ny];   //tau[a][h][l]�����(h,l)������ı��ݣ�taʱ��
static double cellv[2][nx][ny][2];  //cellv[a][h][l][k]Ϊ��(h,l)������ĵ�k���ٶȷ�����taʱ��
static double E[2][nx][ny]; //E[a][h][l]Ϊ(h,l)���������������taʱ��
static double cellP[2][nx][ny]; //cellP[a][h][l]Ϊ(h,l)�����ѹ����taʱ��
static double m[nx][ny];    //��(h,l)��������������غ�
static double e[2][nx][ny]; //��(h,l)����������ܣ�taʱ��
static double c[2][nx][ny]; //��(h,l)����������٣�taʱ��
static double A[2][nx][ny]; //��(h,l)������������taʱ��
static double gamma[nx][ny];    //��(h,l)�������״̬���̲������غ�

//���ڻ��ֵ���

vector<int> neighbornode[n+1];   //��(i,j)���ڵ���������ڽڵ�
vector<int> neighborcell[n+1];   //��(i,j)���������������
vector<double> L[2][n+1];    //��(i,j)���ڵ�ĵ�k���ߵı߳�
vector<double> Nx[2][n+1];   //��(i,j)���ڵ�ĵ�k�����ⷽ��x����
vector<double> Ny[2][n+1];   //��(i,j)���ڵ�ĵ�k�����ⷽ��y����
vector<double> Pstar[n+1][2];    //��(i,j)���ڵ��k�������ߵ�ѹ��
//******************���鶨�����**************//

enum time {tlast, tnext};   //����ö�����ͣ�tlastΪ��һ��ʱ��㣬tnestΪ��һ��ʱ��㣬�ֱ����0��1

//**************��Ҫ�ĺ�������****************//
void initialmesh(); //��ʼ��

void internalsolver(int t_l, int qm);  //�����ڵ���ٶ�
void getinternalPstar(int t_l, int qm);   //�����ڵ���ϵ�ѹ��
void boundarysolverNP(int t_l, int qm);  //����߽�㣨��ƽ�У����ٶ�
void boundarysolverP(int t_l, int qm);    //����߽�㣨ƽ�У����ٶ�
void anglesolver1(int t_l, int qm); //(nx,0)����ٶ�
void anglesolver2(int t_l, int qm); //(nx,ny)����ٶ�
void boundarywithP(int t_l, int qm);    //�߽�������ѹ���������ٶ�
void boundarysolver2(int t_l, int qm);  //����߽������Ϊ2ʱ���ٶ�
void getboundaryPstar(int t_l, int qm);   //����߽���ϱ��ϵ�ѹ��
void nodalsolver(int t_l); //������һʱ������ڵ���ٶȣ����ݽڵ�����ѡ���������������е�һ��
double choosedt(int t_l, double dtn);   //����ʱ�䲽��
void updategeometric(int t_l, int t_n, double dt);  //������ʱ���ļ�����
void updatephysicalvariables(int t_l, int t_n, double dt);   //������ʱ����������
void getstates(int t_n);    //����״̬�������ʣ�±���

//***************�����������****************//

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
    
    fstream fv, frho, fp, fe;   //�ĸ�������

    //�ļ�·��
    const char* filev = "cell-centered_Lagrangian_Scheme\\Noh\\output\\V.txt";    //��¼��Ԫ�ٶ�
    const char* filerho = "cell-centered_Lagrangian_Scheme\\Noh\\output\\rho.txt";    //��¼��Ԫ�ܶ�
    const char* filep = "cell-centered_Lagrangian_Scheme\\Noh\\output\\P.txt";    //��¼��Ԫѹ��
    const char* filee = "cell-centered_Lagrangian_Scheme\\Noh\\output\\e.txt";    //��¼��Ԫ����

    remove(filev);
    remove(filerho);
    remove(filep);
    remove(filee);
    
    do{
        //�ȼ�¼���ݣ����̫���˲����ˣ�600M��

        /*t
        //��¼�ٶ�
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

        //��¼��Ԫ�ܶ�
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

        //��¼��Ԫѹ��
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

        //��¼��������
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
        //��������
        nodalsolver(tl);

        dt = choosedt(tl, dt);
        t = t+ dt;

        updategeometric(tl,tn,dt);
        updatephysicalvariables(tl,tn,dt);
        getstates(tn);

        ttemp = tl;
        tl = tn;
        tn = ttemp; //ԭ��tlʱ�̵���������tnʱ�̵�������һ�μ���ʱtn��Ϊ��һʱ�̵�����
        cout<<t<<"\t"<<dt<<endl<<endl;
    }while(t < 0.6);

    //���濪ʼ��ͼ
    {
        //�ļ���
        const char* plotv = "cell-centered_Lagrangian_Scheme\\Noh\\output\\V.plt";    //��Ԫ�ٶ�
        const char* plotrho = "cell-centered_Lagrangian_Scheme\\Noh\\output\\rho.plt";    //��Ԫ�ܶ�
        const char* plotp = "cell-centered_Lagrangian_Scheme\\Noh\\output\\P.plt";    //��Ԫѹ��
        const char* plote = "cell-centered_Lagrangian_Scheme\\Noh\\output\\e.plt";    //��Ԫ����
        const char* point = "cell-centered_Lagrangian_Scheme\\Noh\\output\\mesh.plt";

        remove(plotv);
        remove(plotrho);
        remove(plotp);
        remove(plote);
        remove(point);

        fstream fpoint;
    
        //�Ȼ��ٶȣ���Ϊ��sod shock������y���ٶ�Ϊ0��ֻ�軭x���ٶ�
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

        //�ٻ��ܶȣ�����sod shock��һά���⣬����ֻҪ��¼x����ͬ
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

        //�ٻ�ѹ��
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

        //�ٻ�����
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

        //�ٻ�����
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

//***********ע�����������ʼ�������ٶȵ�ʱ������²�ͬ***********//
void initialmesh()
{
    int i, j, h, l, q, qt;
    double temp, temp1, temp2;

    //�ȳ�ʼ���ڵ������
    for (q=0; q<=n; q++)
    {
        //�ڵ�����
        if (q == 0)
        {
            i = 0;
            j = 0;
        }
        else{
            i = (q - 1) / (ny + 1) + 1;
            j = (q - 1) % (ny + 1);
        }

            loc[0][q][0] = (i * hx) * cos(j * hy);   //x����
            loc[0][q][1] = (i * hx) * sin(j * hy);   //y����

            //�ж��Ƿ��Ǳ߽��
            if (i == 0)
            {
                boundary[q] = 1; //���㣬��ʱ�߽������ɷ����ٶȸ���
            }
            else if ( i == nx)
            {
                boundary[q] = 2; //�����ϵĵ㣬��ʱbc��ѹ������
            }
            else if (j == 0 || j == ny)
            {
                boundary[q] = 1; //���ϵĵ㣬��ʱbc�ɷ����ٶȸ���
            }
            else{
                boundary[q] = 0;
            }   //�߽���ʼ�����
    }
    //�ڵ��ʼ�����

    //��ʼ��������
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

                e[0][h][l] = cellP[0][h][l] * tau[0][h][l] / (gamma[h][l] - 1);  //����

                E[0][h][l] = 0.5 * (cellv[0][h][l][0] * cellv[0][h][l][0] + cellv[0][h][l][1] * cellv[0][h][l][1]);
                E[0][h][l] = E[0][h][l] + e[0][h][l];   //������
            
                c[0][h][l] = (gamma[h][l] - 1) * e[0][h][l];
                c[0][h][l] = sqrt(c[0][h][l] ); //����

                A[0][h][l] = 0.5 * ( (h+1) * hx ) * ((h+1) * hx) * hy;
                A[0][h][l] = A[0][h][l] - 0.5 * ( h * hx ) * (h * hx) * hy;   //���
                m[h][l] = A[0][h][l] / tau[0][h][l];   //����
            }
            
        }
    }
    //��������ʼ�����

    //��ʼ��������

    //�ȶԶ�����г�ʼ��
    {
        //�ȷ����ڵ�
        for (l = 0; l<=ny; l++)
        {
            neighbornode[0].push_back(l+1);
        }

        //�ٷ������ڵ�Ԫ
        for (l=0; l<ny; l++)
        {
            neighborcell[0].push_back(l);
        }
        neighborcell[0].push_back(-1);

        //�ټ��㳤��
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

        //�ټ����ⷨ����
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

    //�ٶ�ʣ�½ڵ��ʼ��
    int* ntemp = new int [4];
    int* ctemp = new int [4];
    for (q=1; q<=n; q++)
    {

        i = (q - 1) / (ny + 1) + 1;
        j = (q - 1) % (ny + 1);
        
        //����ÿ������ڵ�
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
        
        //�ٰ���ʱ���ڵ�������������
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
        //bΪ�ڵ���������ʱ��һ���ڵ�ı��
        bb = b;

        b = bb;
        //��¼��(i,j)���ڵ���ڵ���
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
        //��¼�ڵ������

        //��¼������
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
            
        //�ж���������
        for (a = 0; a < 4; a++)
        {
            b = a + 1; 
            if (b>3)
            {
                b = b - 4;
            }

            if (ntemp[a] < 0)
            {
                ctemp[a] = -2;  //-2��ʾ����Ҫ����
            }
            else{
                if (ntemp[b] < 0 )
                {
                    ctemp[a] = -1;  //-1��ʾ������
                }
            }
        }

        //��¼����
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
        //��¼�������

        //����ÿ���㵽���ڽڵ�ĳ��Ⱥͷ�����
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
    //��������ʼ�����

    return;
}

void internalsolver(int t_l, int qm)
{
    int h, l, q, k, k2;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22, vk;
    
    //��ʼ������
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


    //������ֵ
    for (k=0; k< neighbornode[qm].size(); k++)
    {
        //�ȼ������ϵ��

        //k2 = k-1��ע����ѭ��ָ��
        {
            k2 = k-1;
            if ( k == 0)
            {
                k2 = neighbornode[qm].size() - 1;
            }
        }

        //ע���ʱ���ڵ�����ÿ������������

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            q = neighborcell[qm][k2];  //(i,j)��ĵ�k�������������

            h = q / ny;
            l = q % ny;

            alpha2 = c[t_l][h][l] / tau[t_l][h][l];
            p2 = cellP[t_l][h][l];
            v21 = cellv[t_l][h][l][0];
            v22 = cellv[t_l][h][l][1];
        }
        
        //����alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            q = neighborcell[qm][k];

            h = q / ny;
            l = q % ny;

            alpha1 = c[t_l][h][l] / tau[t_l][h][l];
            p1 = cellP[t_l][h][l];
            v11 = cellv[t_l][h][l][0];
            v12 = cellv[t_l][h][l][1];
        }

        //����ϵ������A
        {
            A[0][0] = A[0][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Nx[t_l][qm][k];
            A[0][1] = A[0][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][0] = A[1][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][1] = A[1][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Ny[t_l][qm][k] * Ny[t_l][qm][k];
        }

        //����vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][qm][k] + v22 * Ny[t_l][qm][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][qm][k] + v12 * Ny[t_l][qm][k] );
            vk = vk / (alpha1 + alpha2);
        }

        //�����Ҷ˾���B
        {
            B[0] = B[0] + L[t_l][qm][k] * ( alpha1 + alpha2) * vk * Nx[t_l][qm][k];
            B[1] = B[1] + L[t_l][qm][k] * ( alpha1 + alpha2) * vk * Ny[t_l][qm][k];
        }

    }


    //�ⷽ�̣���ά����ֱ���������
    {
        double** Ainv = new double* [2];
        for (h=0; h<2; h++)
        {
            Ainv[h] = new double [2];
        }

        alpha1 = A[0][0] * A[1][1] - A[0][1] * A[1][0]; //alpha1ΪA������ʽ

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
        //k2 = k-1��ע����ѭ��ָ��
        {
            k2 = k-1;
            if ( k == 0)
            {
                k2 = neighbornode[qm].size() - 1;
            }
        }

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            q = neighborcell[qm][k2];  //(i,j)��ĵ�k�������������

            h = q / ny;
            l = q % ny;

            alpha2 = c[t_l][h][l] / tau[t_l][h][l];
            p2 = cellP[t_l][h][l];
            v21 = cellv[t_l][h][l][0];
            v22 = cellv[t_l][h][l][1];
        }
        
        //����alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            q = neighborcell[qm][k];

            h = q / ny;
            l = q % ny;

            alpha1 = c[t_l][h][l] / tau[t_l][h][l];
            p1 = cellP[t_l][h][l];
            v11 = cellv[t_l][h][l][0];
            v12 = cellv[t_l][h][l][1];
        }

        //��P*_{q,k}^{k-1}
        {
            q = neighborcell[qm][k2];
            h = q / ny;
            l = q % ny;

            Pstar[qm][0][k] = 0;

            Pstar[qm][0][k] = Pstar[qm][0][k] + (nodev[qm][0] - cellv[t_l][h][l][0]) * Nx[t_l][qm][k];
            Pstar[qm][0][k] = Pstar[qm][0][k] + (nodev[qm][1] - cellv[t_l][h][l][1]) * Ny[t_l][qm][k];
            Pstar[qm][0][k] = - Pstar[qm][0][k] * alpha2 + p2;
        }

        //��P*_{q,k}^k
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

void boundarysolverNP(int t_l, int qm)    //��ʱ�㷨Ϊ�����߽��ⷨ���ٶȵĽⷨ
{
    int h, l, k;
    double alpha1;
    
    //��ʼ������
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

    //������ֵ
    //����ϵ������A
    {
        k = neighbornode[qm].size() - 1;
        A[0][0] = - Nx[t_l][qm][0];
        A[0][1] = - Ny[t_l][qm][0];
        A[1][0] = Nx[t_l][qm][k]; //  �������㷨��һ��-��
        A[1][1] = Ny[t_l][qm][k]; //�������㷨��һ��-��
    }

    //�����Ҷ˾���B
    {
        B[0] = vn;
        B[1] = vn;
    }

    //�ⷽ�̣���ά����ֱ���������
    {
        double** Ainv = new double* [2];
        for (h=0; h<2; h++)
        {
            Ainv[h] = new double [2];
        }

        alpha1 = A[0][0] * A[1][1] - A[0][1] * A[1][0]; //alpha1ΪA������ʽ

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
    
    //��ʼ������
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

    //������ֵ
    for (k=0; k< neighbornode[qm].size(); k++)
    {
        //�ȼ������ϵ��

        //k2 = k-1��ע��˴�����ѭ��ָ��
        k2 = k - 1;

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            if (  k == 0 )  //-1�ڵ��Ӧ������Ϊ������
            {
                alpha2 = 0;
                p2 = 0;
                v21 = 0;
                v22 = 0;
            }
            else{
                q = neighborcell[qm][k2];  //(i,j)��ĵ�k�������������

                h = q / ny;
                l = q % ny;

                alpha2 = c[t_l][h][l] / tau[t_l][h][l];
                p2 = cellP[t_l][h][l];
                v21 = cellv[t_l][h][l][0];
                v22 = cellv[t_l][h][l][1];
            }
            
        }
        
        //����alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            if ( k == neighbornode[qm].size() - 1 ) //��k���ڵ��Ӧ������Ϊ������
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

        //����ϵ������A
        {
            A[0][0] = A[0][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Nx[t_l][qm][k];
            A[0][1] = A[0][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][0] = A[1][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][1] = A[1][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Ny[t_l][qm][k] * Ny[t_l][qm][k];
            //��ֻ����k�йصĸ�ֵ
        }

        //����vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][qm][k] + v22 * Ny[t_l][qm][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][qm][k] + v12 * Ny[t_l][qm][k] );
        }

        //�����Ҷ˾���B
        {
            B[0] = B[0] + L[t_l][qm][k] * vk * Nx[t_l][qm][k];
            B[1] = B[1] + L[t_l][qm][k] * vk * Ny[t_l][qm][k];
        }

    }
    
    //��ʣ�µĸ�ֵ
    {
        k2 = neighbornode[qm].size() - 1;

        A[0][2] = - L[t_l][qm][0] * Nx[t_l][qm][0] + L[t_l][qm][k2] * Nx[t_l][qm][k2];
        A[2][0] = A[0][2];
        A[1][2] = - L[t_l][qm][0] * Ny[t_l][qm][0] + L[t_l][qm][k2] * Ny[t_l][qm][k2];
        A[2][1] = A[1][2];
        A[2][2] = 0;

        B[2] = L[t_l][qm][0] * vn + L[t_l][qm][k2] * vn;
    }

    //�������Է����飬ע��߽��ⷨ��������ʱ������ά�ģ�ֻȡ�ٶȵĲ���

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

    
    //��A��ֵ
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

    //��B��ֵ
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

    //�ⷽ��
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

    
    //��A��ֵ
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

    //��B��ֵ
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

    //�ⷽ��
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

    //�ȸ�����ֵ
    for (k=0; k< neighbornode[qm].size(); k++)
    {
        //�ȼ������ϵ��

        //k2 = k-1��ע��˴�����ѭ��ָ��
        k2 = k - 1;

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            if (  k == 0 )  //-1�ڵ��Ӧ������Ϊ������
            {
                alpha2 = 0;
                p2 = P_0;
                v21 = 0;
                v22 = 0;
            }
            else{
                q = neighborcell[qm][k2];  //(i,j)��ĵ�k�������������

                h = q / ny;
                l = q % ny;

                alpha2 = c[t_l][h][l] / tau[t_l][h][l];
                p2 = cellP[t_l][h][l];
                v21 = cellv[t_l][h][l][0];
                v22 = cellv[t_l][h][l][1];
            }
            
        }
        
        //����alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            if ( k == neighbornode[qm].size() - 1 ) //��k���ڵ��Ӧ������Ϊ������
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

        //����ϵ������A
        {
            A[0][0] = A[0][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Nx[t_l][qm][k];
            A[0][1] = A[0][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][0] = A[1][0] + L[t_l][qm][k] * ( alpha1 + alpha2) * Nx[t_l][qm][k] * Ny[t_l][qm][k];
            A[1][1] = A[1][1] + L[t_l][qm][k] * ( alpha1 + alpha2) * Ny[t_l][qm][k] * Ny[t_l][qm][k];
            //��ֻ����k�йصĸ�ֵ
        }

        //����vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][qm][k] + v22 * Ny[t_l][qm][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][qm][k] + v12 * Ny[t_l][qm][k] );
        }

        //�����Ҷ˾���B
        {
            B[0] = B[0] + L[t_l][qm][k] * vk * Nx[t_l][qm][k];
            B[1] = B[1] + L[t_l][qm][k] * vk * Ny[t_l][qm][k];
        }

    }

     //�ⷽ�̣���ά����ֱ���������
    {
        double** Ainv = new double* [2];
        for (h=0; h<2; h++)
        {
            Ainv[h] = new double [2];
        }

        alpha1 = A[0][0] * A[1][1] - A[0][1] * A[1][0]; //alpha1ΪA������ʽ

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
    
    if ( qm == (nx - 1)*(ny+1) + 1 )    //qΪ(nx,0)�����
    {
        anglesolver1(t_l, qm);
    }
    else if ( qm == n ) //(nx, ny)�����
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

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            if ( k == 0 )
            {
                alpha2 = 0;
                p2 = 0;
                v21 = 0;
                v22 = 0;
            }
            else{
                q = neighborcell[qm][k2];  //(i,j)��ĵ�k�������������

                h = q / ny;
                l = q % ny;

                alpha2 = c[t_l][h][l] / tau[t_l][h][l];
                p2 = cellP[t_l][h][l];
                v21 = cellv[t_l][h][l][0];
                v22 = cellv[t_l][h][l][1];
            }
        }
        
        //����alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
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

        //��P*_{q,k}^{k-1}
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

        //��P*_{q,k}^k
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

        //�ж��Ƿ��Ǳ߽�㣬���ݲ�ͬ�����
        if (boundary[qm] == 0)    //���ڵ�
        {
            internalsolver(t_l,qm);    //�����(i,j)������ٶ�
            getinternalPstar(t_l,qm);  //�����(i,j)����ı��ϵ�ѹ��
        }
        else if(boundary[qm] == 1)  //�߽������ɷ����ٶȸ���
        {
            //�жϱ߽��Ƿ�ƽ��
            double e = 1e-10;

            q = neighbornode[qm].size()-1;
                
            x1 = Nx[t_l][qm][0];
            y1 = Ny[t_l][qm][0];
            x2 = Nx[t_l][qm][q];
            y2 = Ny[t_l][qm][q];

            temp = x1 * y2 - x2 * y1;
            temp = abs(temp);

            if (temp >e)    //��ƽ��ʱ
            {
                boundarysolverNP(t_l,qm);
            }
            else{
                boundarysolverP(t_l,qm);
            }
            getboundaryPstar(t_l,qm);
        }
        else{   //�߽�����ͬʱ�ɷ����ٶȺ�ѹ������
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
                //�ȼ�¼(h,l)����ÿ������ı��
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

                //ѡ��ÿ����������С�ı߳�
                {
                    lambda = 10000;
                    for (k=0; k<4; k++)
                    {
                        //qmΪ(h,l)�ĵ�k������
                        qm = vertex[k];
                        
                        //�����������㵽(i,j)�����̾���
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

                //ѡȡte
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
                //�ȼ�¼(h,l)����ÿ������ı��
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

                //����ÿ�����������仯��
                {
                    temp = 0;   //temp��¼����仯��
                    for (k=0; k<4; k++)
                    {
                        //qmΪ(h,l)�ĵ�k������
                        qm = vertex[k];

                        //k2 = k + 1, ע����ѭ��ָ��
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

                //ѡȡtv
                
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

    //ѡȡdt
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

    //�ȸ��½ڵ�λ��
    {
        for (qm=0; qm<= n; qm++)
        {
            loc[t_n][qm][0] = loc[t_l][qm][0] + dt * nodev[qm][0];
            loc[t_n][qm][1] = loc[t_l][qm][1] + dt * nodev[qm][1];
        }
    }

    //�ٸ����������
    {
        int* vertex = new int [4];

        for (h=0; h<nx; h++)
        {
            for (l=0; l <ny; l++)
            {
                A[t_n][h][l] = 0;
                //�ȼ�¼(h,l)����ÿ������ı��
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
                    //k2 = k+1,ע����ѭ��ָ��
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

    //�ٸ�������ı߳����ⷨ����
    {
        
        for (qm=0; qm <= n; qm++)
        {
            for (k=0; k<neighbornode[qm].size(); k++)
            {
                qm2 = neighbornode[qm][k];

                //�ȸ��±߳�
                temp = ( loc[t_n][qm2][0] - loc[t_n][qm][0] ) * ( loc[t_n][qm2][0] - loc[t_n][qm][0] );
                temp = temp + ( loc[t_n][qm2][1] - loc[t_n][qm][1] ) * ( loc[t_n][qm2][1] - loc[t_n][qm][1] );
                temp = sqrt(temp);

                L[t_n][qm][k] = temp;

                //�ٸ����ⷨ����
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
            //��¼�������
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

            //�ȸ��±���
            {
                temp = 0;
                for (k=0; k< 4; k++)
                {
                    qm = vertex[k];   //qmΪ(h,l)�ĵ�k������

                    k2 = k + 1;
                    if (k2 > 3)
                    {
                        k2 = k2 - 4;
                    }

                    qm2 = vertex[k2];  //qm2Ϊ(h,l)�ĵ�k+1������

                    //ȷ�ϵ�k���߳��ͷ�����
                    for (r=0; r<neighbornode[qm].size(); r++)
                    {
                        if ( neighbornode[qm][r] == vertex[k2] )  //ȷ��k+1������k����ĵڼ����ڵ�
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
                
                temp = - temp * dt * 0.25;  //�����������ⷨ������N�ķ����෴�õ�
                temp = temp / m[h][l];
                
                temp = temp + tau[t_l][h][l];

                tau[t_n][h][l] = temp;
            }

            //�ٸ����ٶ�
            {
                temp1 = 0;
                temp2 = 0;

                for (k=0; k< 4; k++)
                {
                    qm = vertex[k];

                    //����L_{r,r+1}^n * P_{r,r+1/2}^{*,i} * N_{r,r+1}^n
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
                                //ע�������ⷨ������N�ķ����ϵ
                                temp1 = temp1 - L[t_l][qm][r] * Pstar[qm][1][r] * Nx[t_l][qm][r];
                                temp2 = temp2 - L[t_l][qm][r] * Pstar[qm][1][r] * Ny[t_l][qm][r];
                                break;
                            }
                        }
                    }
                    
                    //����L_{r,r+1}^n * P_{r+1/2, r+1}^{*,i} * N_{r,r+1}^n
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
            
            //�ٸ�������
            {
                temp = 0;

                for (k=0; k<4; k++)
                {
                    qm = vertex[k];

                    //�ȼ���L_{r,r+1}^n * P_{r,r+1/2}^{*,i} * V*_r * N_{r,r+1}^n
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
                                temp = temp - L[t_l][qm][r] * Pstar[qm][1][r] * temp1;  //ע����������ⷨ������N�ķ����ϵ
                                break;
                            }
                        }

                    }

                    //����L_{r,r+1}^n * P_{r+1/2, r+1}^{*,i} * V*_{r+1} * N_{r,r+1}^n
                    //ע��ѭ��ָ��ı仯
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
            //����ÿ����Ԫ������
            {
                temp = cellv[t_n][h][l][0] * cellv[t_n][h][l][0] + cellv[t_n][h][l][1] * cellv[t_n][h][l][1];
                e[t_n][h][l] = E[t_n][h][l] - 0.5 * temp;
            }

            //�ٸ���״̬������Ԫѹ��
            {
                cellP[t_n][h][l] = ( gamma[h][l] - 1 ) * e[t_n][h][l];
                cellP[t_n][h][l] = cellP[t_n][h][l] / tau[t_n][h][l]; 
            }

            //�������
            {
                c[t_n][h][l] = (gamma[h][l] - 1) * e[t_n][h][l];
                c[t_n][h][l] = sqrt(c[t_n][h][l]);
            }
        }
    }

    return ;
}
