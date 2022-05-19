#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "initial.h"

using namespace std;

//************���涨����Ҫ�õ�������***********//

//���ڽڵ����

static double loc[2][nx+1][ny+1][2]; //loc[a][i][j][k]Ϊ��(i,j)���ڵ�taʱ�̵ĵ�k��������ta=tn��tn+1
static double nodev[nx+1][ny+1][2];    //nodev[i][j][k]Ϊ��(i,j)����ĵ�k���ٶȷ����������㷨ֻ��Ҫһ������洢
static bool boundary[nx+1][ny+1];   //����0�������ڵ㣬�����Ǳ߽��

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

vector<int> neighbornode[nx+1][ny+1];   //��(i,j)���ڵ���������ڽڵ�
vector<int> neighborcell[nx+1][ny+1];   //��(i,j)���������������
vector<double> L[2][nx+1][ny+1];    //��(i,j)���ڵ�ĵ�k���ߵı߳�
vector<double> Nx[2][nx+1][ny+1];   //��(i,j)���ڵ�ĵ�k�����ⷽ��x����
vector<double> Ny[2][nx+1][ny+1];   //��(i,j)���ڵ�ĵ�k�����ⷽ��y����
vector<double> Pstar[nx+1][ny+1][2];    //��(i,j)���ڵ��k�������ߵ�ѹ��
//******************���鶨�����**************//

enum time {tlast, tnext};   //����ö�����ͣ�tlastΪ��һ��ʱ��㣬tnestΪ��һ��ʱ��㣬�ֱ����0��1

//**************��Ҫ�ĺ�������****************//
void initialmesh(); //��ʼ��

void internalsolver(int t_l, int i, int j);  //�����ڵ���ٶ�
void getinternalPstar(int t_l, int i, int j);   //�����ڵ���ϵ�ѹ��
void boundarysolverNP(int t_l, int i, int j);  //����߽�㣨��ƽ�У����ٶ�
void boundarysolverP(int t_l, int i, int j);    //����߽�㣨ƽ�У����ٶ�
void getboundaryPstar(int t_l, int i, int j);   //����߽���ϱ��ϵ�ѹ��
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

    /*for (i=0; i<=nx; i++)
    {
        for (j=0; j<=ny; j++)
        {
            q = i * (ny + 1) + j;
            cout<<q<<"(";
            
            for (h=0; h<L[tn][i][j].size(); h++)
            {
                cout<<L[tn][i][j][h]<<"\t";
            }
            
            if (i < nx && j < ny)
            {
                cout<<cellP[tn][i][j];
            }
            
            cout<<")";
            cout<<"\t";
        }
        cout<<endl;
    }
    */

    fstream fv, frho, fp, fe;   //�ĸ�������

    //�ļ�·��
    const char* filev = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\V.txt";    //��¼��Ԫ�ٶ�
    const char* filerho = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\rho.txt";    //��¼��Ԫ�ܶ�
    const char* filep = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\P.txt";    //��¼��Ԫѹ��
    const char* filee = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\e.txt";    //��¼��Ԫ����

    remove(filev);
    remove(filerho);
    remove(filep);
    remove(filee);
    
    do{

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
        cout<<t<<"\t"<<dt<<endl;
    }while(t < 0.2);

    //���濪ʼ��ͼ
    {
        //�ļ���
        const char* plotv = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\V.plt";    //��Ԫ�ٶ�
        const char* plotrho = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\rho.plt";    //��Ԫ�ܶ�
        const char* plotp = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\P.plt";    //��Ԫѹ��
        const char* plote = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\e.plt";    //��Ԫ����
        const char* point = "cell-centered_Lagrangian_Scheme\\sod_shock\\output\\mesh.plt";
        remove(plotv);
        remove(plotrho);
        remove(plotp);
        remove(plote);
        remove(point);

        fstream fpoint;
    
        //�Ȼ��ٶȣ���Ϊ��sod shock������y���ٶ�Ϊ0��ֻ�軭x���ٶ�
        {
            fv.open(plotv, ios :: out | ios :: app);
            fv<<"VARIABLES="<<"X"<<","<<"V, V_exact"<<endl;
            for (h=0; h<nx; h++)
            {
                double vexact;
                if (loc[tl][h][0][0]<0.2*-1.29099 + 0.5){
                    vexact = 0;
                }
                else if (loc[tl][h][0][0]<0.2*-0.209305 + 0.5)
                {
                    vexact = 0.901408 * (loc[tl][h][0][0] - (0.2*-1.29099 + 0.5)) / ((0.2*-0.209305 + 0.5) -(0.2*-1.29099 + 0.5));
                }
                else if (loc[tl][h][0][0]<0.2*1.90265 + 0.5)
                {
                    vexact = 0.901408;
                }
                else
                {
                    vexact = 0;
                }
                fv<<"\t"<<loc[tl][h][0][0]<<"\t"<<cellv[tl][h][0][0]<<"\t"<<vexact<<endl;
            }
            fv.close();
        }

        //�ٻ��ܶȣ�����sod shock��һά���⣬����ֻҪ��¼x����ͬ
        {
            frho.open(plotrho, ios :: out | ios :: app);
            frho<<"VARIABLES="<<"X"<<","<<"density, d_exact"<<endl;
            for (h=0; h<nx; h++)
            {
                double rexact;
                if (loc[tl][h][0][0]<0.2*-1.29099 + 0.5){
                    rexact = rho_l;
                }
                else if (loc[tl][h][0][0]<0.2*-0.209305 + 0.5)
                {
                    rexact = (0.424722-rho_l) * (loc[tl][h][0][0] - (0.2*-1.29099 + 0.5)) / ((0.2*-0.209305 + 0.5) -(0.2*-1.29099 + 0.5));
                    rexact = rexact + rho_l;
                }
                else if (loc[tl][h][0][0]<0.2*0.901408 + 0.5)
                {
                    rexact = 0.424722;
                }
                else if (loc[tl][h][0][0]<0.2*1.90265 + 0.5)
                {
                    rexact = 0.237536;
                }
                else
                {
                    rexact = rho_r;
                }
                temp = 1 / tau[tl][h][0];
                frho<<"\t"<<loc[tl][h][0][0]<<"\t"<<temp<<"\t"<<rexact<<endl;
            }
            frho.close();
        }

        //�ٻ�ѹ��
        {
            fp.open(plotp, ios :: out | ios :: app);
            fp<<"VARIABLES="<<"X"<<","<<"P, P_exact"<<endl;
            for (h=0; h< nx; h++)
            {
                double pexact;
                if (loc[tl][h][0][0]<0.2*-1.29099 + 0.5){
                    pexact = P_l;
                }
                else if (loc[tl][h][0][0]<0.2*-0.209305 + 0.5)
                {
                    pexact = (0.314383- P_l) * (loc[tl][h][0][0] - (0.2*-1.29099 + 0.5)) / ((0.2*-0.209305 + 0.5) -(0.2*-1.29099 + 0.5));
                    pexact = pexact + P_l;
                }
                else if (loc[tl][h][0][0]<0.2*1.90265 + 0.5)
                {
                    pexact = 0.314383;
                }
                else
                {
                    pexact = P_r;
                }
                fp<<"\t"<<loc[tl][h][0][0]<<"\t"<<cellP[tl][h][0]<<"\t"<<pexact<<endl;
            }
            fp.close();
        }

        //�ٻ�����
        {
            fe.open(plote, ios :: out | ios :: app);
            fe<<"VARIABLES="<<"X"<<","<<"e, e_exact"<<endl;
            for (h=0; h< nx; h++)
            {
                double eexact;
                if (loc[tl][h][0][0]<0.2*-1.29099 + 0.5){
                    eexact = P_l/((gamma_l-1)*rho_l);
                }
                else if (loc[tl][h][0][0]<0.2*-0.209305 + 0.5)
                {
                    double el, er;
                    el = P_l/((gamma_l-1)*rho_l);
                    er = 0.314383 / ((gamma_l-1)*0.424722);
                    eexact = (er-el) * (loc[tl][h][0][0]-(0.2*-1.29099 + 0.5)) / (0.2*-0.209305 + 0.5 - (0.2*-1.29099 + 0.5));
                    eexact = eexact + el;
                }
                else if (loc[tl][h][0][0]<0.2*0.901408 + 0.5)
                {
                    double pl, rhl;
                    pl = 0.314383;
                    rhl = 0.424722;
                    eexact = pl / ((gamma_l-1)*rhl);
                }
                else if (loc[tl][h][0][0]<0.2*1.90265 + 0.5)
                {
                    double pl, rhl;
                    pl = 0.314383;
                    rhl = 0.237536;
                    eexact = pl / ((gamma_r-1)*rhl);
                }
                else
                {
                    eexact = P_r/((gamma_r-1)*rho_r);
                }
                fe<<"\t"<<loc[tl][h][0][0]<<"\t"<<e[tl][h][0]<<"\t"<<eexact<<endl;
            }
            fe.close();
        }
        
        //�ٻ�����
        {
            fpoint.open(point, ios :: out | ios :: app);
            fpoint<<"VARIABLES="<<"X"<<","<<"Y"<<endl;
            fpoint<<"ZONE N="<<n<<", E="<<nx*ny<<", F=FEPOINT, ET=QUADRILATERAL"<<endl;
            for (i=0; i<=nx; i++)
            {
                for (j=0; j<=ny; j++)
                {
                    fpoint<<"\t"<<loc[tl][i][j][0]<<"\t"<<loc[tl][i][j][1]<<endl;
                }
            }
                double ctemp[4];
                
            for (h=0; h<nx; h++)
            {
                for (l=0; l<ny; l++)
                {
                    ctemp[0] = h * (ny+1) + l;
                    ctemp[1] = h * (ny+1) + l + 1;
                    ctemp[2] = (h+1) * (ny+1) + l+1;
                    ctemp[3] = (h+1) * (ny+1) + l;
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

void initialmesh()
{
    //�ȳ�ʼ���ڵ������
    int i, j, h, l, q;
    double temp, temp1, temp2;

    for (i=0; i<= nx; i++)
    {
        for (j=0; j<=ny; j++)
        {
            loc[0][i][j][0] = hx * i;   //x����
            loc[0][i][j][1] = hy * j;   //y����

            //�ж��Ƿ��Ǳ߽��
            if (i == 0 || i == nx)
            {
                boundary[i][j] = 1;
            }
            else if (j == 0 || j == ny)
            {
                boundary[i][j] = 1;
            }
            else{
                boundary[i][j] = 0;
            }   //�߽���ʼ�����
        }
    }
    //�ڵ��ʼ�����

    //��ʼ��������
    for (h=0; h<nx; h++)
    {
        for (l=0; l<ny; l++)
        {
            //������ֵ
            if ( (h*hx) < x/2 )
            {
                tau[0][h][l] = 1/ rho_l;    //����
                cellv[0][h][l][0] = u_l;    //��Ԫ�ٶ�
                cellP[0][h][l] = P_l;   //��Ԫѹ��
                gamma[h][l] = gamma_l;  //��Ԫ״̬���̲���
            }
            else{
                tau[0][h][l] = 1/ rho_r;
                cellv[0][h][l][1] = u_r;
                cellP[0][h][l] = P_r;
                gamma[h][l] = gamma_r;
            }
            
            {
                e[0][h][l] = cellP[0][h][l] * tau[0][h][l] / (gamma[h][l] - 1);  //����

                E[0][h][l] = 0.5 * (cellv[0][h][l][0] * cellv[0][h][l][0] + cellv[0][h][l][1] * cellv[0][h][l][1]);
                E[0][h][l] = E[0][h][l] + e[0][h][l];   //������
            
                c[0][h][l] = (gamma[h][l] - 1) * e[0][h][l];
                c[0][h][l] = sqrt(c[0][h][l] ); //����

                A[0][h][l] = hx * hy;   //���
                m[h][l] = A[0][h][l] / tau[0][h][l];   //����
            }
            
        }
    }
    //��������ʼ�����

    //��ʼ��������
    int* ntemp = new int [4];
    int* ctemp = new int [4];
    for (i=0; i<=nx; i++)
    {
        for (j=0; j<=ny; j++)
        {
            //����ÿ������ڵ�
            {
                if (j-1 < 0)    //(i,j-1)
                {
                    ntemp[0] = -1;
                }
                else{
                    ntemp[0] = i * (ny + 1) + j - 1;
                }
                if (i+1 > nx)   //(i+1,j)
                {
                    ntemp[1] = -1;
                }
                else{
                    ntemp[1] = (i+1) * (ny + 1 ) + j;
                }
                if (j+1 > ny)   //(i,j+1)
                {
                    ntemp[2] = -1;
                }
                else{
                    ntemp[2] = i * (ny + 1) + j + 1;
                }
                if (i-1 < 0)    //(i-1,j)
                {
                    ntemp[3] = -1;
                }
                else{
                    ntemp[3] = (i-1) * (ny + 1) + j;
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
                    neighbornode[i][j].push_back(ntemp[b]);
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
                ctemp[2] = (i-1) * ny + j;
                ctemp[3] = (i-1) * ny + j - 1;
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
                    neighborcell[i][j].push_back(ctemp[b]);
                }
                b++;
                if (b>3)
                {
                    b = b-4;
                }
            }
            //��¼�������

            //����ÿ���㵽���ڽڵ�ĳ��Ⱥͷ�����
            for (bb = 0; bb < neighbornode[i][j].size(); bb++)
            {
                q = neighbornode[i][j][bb];
                h = q / (ny + 1);
                l = q % (ny + 1); //(h,l)Ϊ�ڵ�����

                temp1 = loc[0][h][l][0] - loc[0][i][j][0];
                temp2 = loc[0][h][l][1] - loc[0][i][j][1];
                temp = temp1 * temp1 + temp2 * temp2;
                temp = sqrt(temp);
                
                L[0][i][j].push_back(temp);
                L[1][i][j].push_back(0);

                Nx[0][i][j].push_back(- temp2 / temp);
                Ny[0][i][j].push_back(temp1 / temp);
                Nx[1][i][j].push_back(0);
                Ny[1][i][j].push_back(0);

                Pstar[i][j][0].push_back(0);
                Pstar[i][j][1].push_back(0);
            }
        }
    }
    delete[] ntemp;
    delete[] ctemp;
    //��������ʼ�����

    return;
}

void internalsolver(int t_l, int i, int j)
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
    for (k=0; k< neighbornode[i][j].size(); k++)
    {
        //�ȼ������ϵ��

        //k2 = k-1��ע����ѭ��ָ��
        {
            k2 = k-1;
            if ( k == 0)
            {
                k2 = neighbornode[i][j].size() - 1;
            }
        }

        //ע���ʱ���ڵ�����ÿ������������

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            q = neighborcell[i][j][k2];  //(i,j)��ĵ�k�������������

            h = q / ny;
            l = q % ny;

            alpha2 = c[t_l][h][l] / tau[t_l][h][l];
            p2 = cellP[t_l][h][l];
            v21 = cellv[t_l][h][l][0];
            v22 = cellv[t_l][h][l][1];
        }
        
        //����alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            q = neighborcell[i][j][k];

            h = q / ny;
            l = q % ny;

            alpha1 = c[t_l][h][l] / tau[t_l][h][l];
            p1 = cellP[t_l][h][l];
            v11 = cellv[t_l][h][l][0];
            v12 = cellv[t_l][h][l][1];
        }

        //����ϵ������A
        {
            A[0][0] = A[0][0] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Nx[t_l][i][j][k] * Nx[t_l][i][j][k];
            A[0][1] = A[0][1] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Nx[t_l][i][j][k] * Ny[t_l][i][j][k];
            A[1][0] = A[1][0] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Nx[t_l][i][j][k] * Ny[t_l][i][j][k];
            A[1][1] = A[1][1] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Ny[t_l][i][j][k] * Ny[t_l][i][j][k];
        }

        //����vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][i][j][k] + v22 * Ny[t_l][i][j][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][i][j][k] + v12 * Ny[t_l][i][j][k] );
            vk = vk / (alpha1 + alpha2);
        }

        //�����Ҷ˾���B
        {
            B[0] = B[0] + L[t_l][i][j][k] * ( alpha1 + alpha2) * vk * Nx[t_l][i][j][k];
            B[1] = B[1] + L[t_l][i][j][k] * ( alpha1 + alpha2) * vk * Ny[t_l][i][j][k];
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
            nodev[i][j][h] = 0;
            for (l=0; l<2; l++)
            {
                nodev[i][j][h] = nodev[i][j][h] + Ainv[h][l] * B[l];
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

void getinternalPstar(int t_l, int i, int j)
{
    int h, l, q, k, k2;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22;

    for (k=0; k<neighbornode[i][j].size(); k++)
    {
        //k2 = k-1��ע����ѭ��ָ��
        {
            k2 = k-1;
            if ( k == 0)
            {
                k2 = neighbornode[i][j].size() - 1;
            }
        }

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            q = neighborcell[i][j][k2];  //(i,j)��ĵ�k�������������

            h = q / ny;
            l = q % ny;

            alpha2 = c[t_l][h][l] / tau[t_l][h][l];
            p2 = cellP[t_l][h][l];
            v21 = cellv[t_l][h][l][0];
            v22 = cellv[t_l][h][l][1];
        }
        
        //����alpha1 = alpha_k, p1 = cellP_k, v1 = cellv_k
        {
            q = neighborcell[i][j][k];

            h = q / ny;
            l = q % ny;

            alpha1 = c[t_l][h][l] / tau[t_l][h][l];
            p1 = cellP[t_l][h][l];
            v11 = cellv[t_l][h][l][0];
            v12 = cellv[t_l][h][l][1];
        }

        //��P*_{q,k}^{k-1}
        {
            q = neighborcell[i][j][k2];
            h = q / ny;
            l = q % ny;

            Pstar[i][j][0][k] = 0;

            Pstar[i][j][0][k] = Pstar[i][j][0][k] + (nodev[i][j][0] - cellv[t_l][h][l][0]) * Nx[t_l][i][j][k];
            Pstar[i][j][0][k] = Pstar[i][j][0][k] + (nodev[i][j][1] - cellv[t_l][h][l][1]) * Ny[t_l][i][j][k];
            Pstar[i][j][0][k] = - Pstar[i][j][0][k] * alpha2 + p2;
        }

        //��P*_{q,k}^k
        {
            q = neighborcell[i][j][k];
            h = q / ny;
            l = q % ny;

            Pstar[i][j][1][k] = 0;

            Pstar[i][j][1][k] = Pstar[i][j][1][k] + (nodev[i][j][0] - cellv[t_l][h][l][0]) * Nx[t_l][i][j][k];
            Pstar[i][j][1][k] = Pstar[i][j][1][k] + (nodev[i][j][1] - cellv[t_l][h][l][1]) * Ny[t_l][i][j][k];
            Pstar[i][j][1][k] = Pstar[i][j][1][k] * alpha1 + p1;
        }
    }

    return;
}

void boundarysolverNP(int t_l, int i, int j)    //��ʱ�㷨Ϊ�����߽��ⷨ���ٶȵĽⷨ
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
        k = neighbornode[i][j].size() - 1;
        A[0][0] = - Nx[t_l][i][j][0];
        A[0][1] = - Ny[t_l][i][j][0];
        A[1][0] = Nx[t_l][i][j][k]; //  �������㷨��һ��-��
        A[1][1] = Ny[t_l][i][j][k]; //�������㷨��һ��-��
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
            nodev[i][j][h] = 0;
            for (l=0; l<2; l++)
            {
                nodev[i][j][h] = nodev[i][j][h] + Ainv[h][l] * B[l];
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

void boundarysolverP(int t_l, int i, int j)
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
    for (k=0; k< neighbornode[i][j].size(); k++)
    {
        //�ȼ������ϵ��

        //k2 = k-1��ע��˴�����ѭ��ָ��
        k2 = k - 1;

        //����alpha2 = alpha_{k-1}, p2 = cellP_{k-1}, v2 = cellv_{k-1}
        {
            if (  k == 0 )
            {
                alpha2 = 0;
                p2 = 0;
                v21 = 0;
                v22 = 0;
            }
            else{
                q = neighborcell[i][j][k2];  //(i,j)��ĵ�k�������������

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
            if ( k == neighbornode[i][j].size() - 1 )
            {
                alpha1 = 0;
                p1 = 0;
                v11 = 0;
                v12 = 0;
            }
            else{
                q = neighborcell[i][j][k];

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
            A[0][0] = A[0][0] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Nx[t_l][i][j][k] * Nx[t_l][i][j][k];
            A[0][1] = A[0][1] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Nx[t_l][i][j][k] * Ny[t_l][i][j][k];
            A[1][0] = A[1][0] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Nx[t_l][i][j][k] * Ny[t_l][i][j][k];
            A[1][1] = A[1][1] + L[t_l][i][j][k] * ( alpha1 + alpha2) * Ny[t_l][i][j][k] * Ny[t_l][i][j][k];
            //��ֻ����k�йصĸ�ֵ
        }

        //����vk
        {
            vk = p2 - p1;
            vk = vk + alpha2 * ( v21 * Nx[t_l][i][j][k] + v22 * Ny[t_l][i][j][k] );
            vk = vk + alpha1 * ( v11 * Nx[t_l][i][j][k] + v12 * Ny[t_l][i][j][k] );
        }

        //�����Ҷ˾���B
        {
            B[0] = B[0] + L[t_l][i][j][k] * vk * Nx[t_l][i][j][k];
            B[1] = B[1] + L[t_l][i][j][k] * vk * Ny[t_l][i][j][k];
        }

    }
    
    //��ʣ�µĸ�ֵ
    {
        k2 = neighbornode[i][j].size() - 1;

        A[0][2] = - L[t_l][i][j][0] * Nx[t_l][i][j][0] + L[t_l][i][j][k2] * Nx[t_l][i][j][k2];
        A[2][0] = A[0][2];
        A[1][2] = - L[t_l][i][j][0] * Ny[t_l][i][j][0] + L[t_l][i][j][k2] * Ny[t_l][i][j][k2];
        A[2][1] = A[1][2];
        A[2][2] = 0;

        B[2] = L[t_l][i][j][0] * vn + L[t_l][i][j][k2] * vn;
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
    
    nodev[i][j][0] = B[0];
    nodev[i][j][1] = B[1];

    for (h=0; h<3; h++)
    {
        delete[] A[h];
    }
    delete[] A;
    delete[] B;
}

void getboundaryPstar(int t_l, int i, int j)
{
    int h, l, q, k, k2;
    double alpha1, alpha2, p1, p2, v11, v12, v21, v22;

    for (k=0; k<neighbornode[i][j].size(); k++)
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
                q = neighborcell[i][j][k2];  //(i,j)��ĵ�k�������������

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
            if ( k == neighbornode[i][j].size() - 1 )
            {
                alpha1 = 0;
                p1 = 0;
                v11 = 0;
                v12 = 0;
            }
            else{
                q = neighborcell[i][j][k];

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
                Pstar[i][j][0][k] = 0;
            }
            else{
                q = neighborcell[i][j][k2];
                h = q / ny;
                l = q % ny;

                Pstar[i][j][0][k] = 0;

                Pstar[i][j][0][k] = Pstar[i][j][0][k] + (nodev[i][j][0] - cellv[t_l][h][l][0]) * Nx[t_l][i][j][k];
                Pstar[i][j][0][k] = Pstar[i][j][0][k] + (nodev[i][j][1] - cellv[t_l][h][l][1]) * Ny[t_l][i][j][k];
                Pstar[i][j][0][k] = - Pstar[i][j][0][k] * alpha2 + p2;
            }
        }

        //��P*_{q,k}^k
        {
            if ( k == neighbornode[i][j].size() - 1 )
            {
                Pstar[i][j][1][k] = 0;
            }
            else{
                q = neighborcell[i][j][k];
                h = q / ny;
                l = q % ny;

                Pstar[i][j][1][k] = 0;

                Pstar[i][j][1][k] = Pstar[i][j][1][k] + (nodev[i][j][0] - cellv[t_l][h][l][0]) * Nx[t_l][i][j][k];
                Pstar[i][j][1][k] = Pstar[i][j][1][k] + (nodev[i][j][1] - cellv[t_l][h][l][1]) * Ny[t_l][i][j][k];
                Pstar[i][j][1][k] = Pstar[i][j][1][k] * alpha1 + p1;
            }
        }
    }

    return;
}

void nodalsolver(int t_l)
{
    int i, j, h, l, q;
    double x1, x2, y1, y2, temp;

    for (i=0; i<=nx; i++)
    {
        for (j=0; j<=ny; j++)
        {
            //�ж��Ƿ��Ǳ߽�㣬���ݲ�ͬ�����
            if (boundary[i][j] == 0)    //���ڵ�
            {
                internalsolver(t_l,i,j);    //�����(i,j)������ٶ�
                getinternalPstar(t_l,i,j);  //�����(i,j)����ı��ϵ�ѹ��
            }
            else{
                //�жϱ߽��Ƿ�ƽ��
                double e = 1e-10;

                q = neighbornode[i][j].size()-1;
                
                x1 = Nx[t_l][i][j][0];
                y1 = Ny[t_l][i][j][0];
                x2 = Nx[t_l][i][j][q];
                y2 = Ny[t_l][i][j][q];

                temp = x1 * y2 - x2 * y1;
                temp = abs(temp);

                if (temp >e)    //��ƽ��ʱ
                {
                    boundarysolverNP(t_l,i,j);
                }
                else{
                    boundarysolverP(t_l,i,j);
                }
                getboundaryPstar(t_l,i,j);
            }

        }
    }
    return;
}

double choosedt(int t_l, double dtn)
{
    double dt, te, tv, lambda, temp;
    int i, j, h, l, k, k2, p, q;

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
                    vertex[0] = h * (ny + 1) + l;
                    vertex[1] = (h + 1) * (ny + 1) + l;
                    vertex[2] = (h + 1) * (ny + 1) + l + 1;
                    vertex[3] = h * (ny + 1) + l + 1;
                }

                //ѡ��ÿ����������С�ı߳�
                {
                    lambda = 10000;
                    for (k=0; k<4; k++)
                    {
                        //(i,j)Ϊ(h,l)�ĵ�k������
                        i = vertex[k] / (ny + 1);
                        j = vertex[k] % (ny + 1);

                        //�����������㵽(i,j)�����̾���
                        k2 = k + 1;
                        if ( k2 > 3)
                        {
                            k2 = k2 - 4;
                        }
                        do{

                            p = vertex[k2] / (ny + 1);
                            q = vertex[k2] % (ny + 1);

                            temp = (loc[t_l][i][j][0] - loc[t_l][p][q][0]) * (loc[t_l][i][j][0] - loc[t_l][p][q][0]);
                            temp = temp + (loc[t_l][i][j][1] - loc[t_l][p][q][1]) * (loc[t_l][i][j][1] - loc[t_l][p][q][1]);
                            temp = sqrt(temp);

                            if (temp < lambda)
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
                    vertex[0] = h * (ny + 1) + l;
                    vertex[1] = (h + 1) * (ny + 1) + l;
                    vertex[2] = (h + 1) * (ny + 1) + l + 1;
                    vertex[3] = h * (ny + 1) + l + 1;
                }

                //����ÿ�����������仯��
                {
                    temp = 0;   //temp��¼����仯��
                    for (k=0; k<4; k++)
                    {
                        //(i,j)Ϊ(h,l)�ĵ�k������
                        i = vertex[k] / (ny + 1);
                        j = vertex[k] % (ny + 1);

                        //k2 = k + 1, ע����ѭ��ָ��
                        k2 = k + 1;
                        if ( k2 > 3)
                        {
                            k2 = k2 - 4;
                        }

                        p = vertex[k2] / (ny + 1);
                        q = vertex[k2] % (ny + 1);

                        lambda = nodev[i][j][0] * loc[t_l][p][q][1];
                        lambda = lambda + nodev[p][q][1] * loc[t_l][i][j][0];
                        lambda = lambda - nodev[p][q][0] * loc[t_l][i][j][1];
                        lambda = lambda - nodev[i][j][1] * loc[t_l][p][q][0];

                        temp = temp + lambda;
                        
                    }
                    temp = temp * 0.5;
                    temp = abs(temp);
                }

                //ѡȡtv
                
                temp = A[t_l][h][l] / temp;

                if (temp < tv)
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
    int i, j, h, l, k, k2, p, q;
    double temp;

    //�ȸ��½ڵ�λ��
    {
        for (i=0; i<= nx; i++)
        {
            for (j=0; j <= ny; j++)
            {
                loc[t_n][i][j][0] = loc[t_l][i][j][0] + dt * nodev[i][j][0];
                loc[t_n][i][j][1] = loc[t_l][i][j][1] + dt * nodev[i][j][1];
            }
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
                    vertex[0] = h * (ny + 1) + l;
                    vertex[1] = (h + 1) * (ny + 1) + l;
                    vertex[2] = (h + 1) * (ny + 1) + l + 1;
                    vertex[3] = h * (ny + 1) + l + 1;
                }
                
                for (k=0; k < 4; k++)
                {
                    //k2 = k+1,ע����ѭ��ָ��
                    k2 = k + 1;
                    if (k2 > 3)
                    {
                        k2 = k2 - 4;
                    }

                    i = vertex[k] / (ny + 1);
                    j = vertex[k] % (ny + 1);

                    p = vertex[k2] / (ny + 1);
                    q = vertex[k2] % (ny + 1);

                    temp = ( loc[t_l][i][j][0] + dt * nodev[i][j][0] ) * ( loc[t_l][p][q][1] + dt * nodev[p][q][1] );
                    temp = temp - ( loc[t_l][i][j][1] + dt * nodev[i][j][1] ) * ( loc[t_l][p][q][0] + dt * nodev[p][q][0] );

                    A[t_n][h][l] = A[t_n][h][l] + temp * 0.5;
                }
            }
        }

        delete[] vertex;
    }

    //�ٸ�������ı߳����ⷨ����
    {
        
        for (i=0; i<= nx; i++)
        {
            for (j=0; j<=ny; j++)
            {
                
                for (k=0; k<neighbornode[i][j].size(); k++)
                {
                    p = neighbornode[i][j][k] / (ny + 1);
                    q = neighbornode[i][j][k] % (ny + 1);

                    //�ȸ��±߳�
                    temp = ( loc[t_n][p][q][0] - loc[t_n][i][j][0] ) * ( loc[t_n][p][q][0] - loc[t_n][i][j][0] );
                    temp = temp + ( loc[t_n][p][q][1] - loc[t_n][i][j][1] ) * ( loc[t_n][p][q][1] - loc[t_n][i][j][1] );
                    temp = sqrt(temp);

                    L[t_n][i][j][k] = temp;

                    //�ٸ����ⷨ����

                    Nx[t_n][i][j][k] = - (loc[t_n][p][q][1] - loc[t_n][i][j][1]) / temp;
                    Ny[t_n][i][j][k] = (loc[t_n][p][q][0] - loc[t_n][i][j][0]) / temp;
                }

            }
        }

    }
    return ;
}

void updatephysicalvariables(int t_l, int t_n, double dt)
{
    int i, j, h, l, k, k2, p, q, r;
    double temp, temp1, temp2;

    int* vertex = new int [4];

    for (h=0; h < nx; h++)
    {
        for (l=0; l<ny; l++)
        {
            //��¼�������
            {
                vertex[0] = h * (ny + 1) + l;
                vertex[1] = (h + 1) * (ny + 1) + l;
                vertex[2] = (h + 1) * (ny + 1) + l + 1;
                vertex[3] = h * (ny + 1) + l + 1;
            }

            //�ȸ��±���
            {
                temp = 0;
                for (k=0; k< 4; k++)
                {
                    i = vertex[k] / (ny + 1);
                    j = vertex[k] % (ny + 1);   //(i,j)Ϊ(h,l)�ĵ�k������

                    k2 = k + 1;
                    if (k2 > 3)
                    {
                        k2 = k2 - 4;
                    }

                    p = vertex[k2] / (ny + 1);
                    q = vertex[k2] % (ny + 1);  //(p,q)Ϊ(h,l)�ĵ�k+1������

                    //ȷ�ϵ�k���߳��ͷ�����
                    for (r=0; r<neighbornode[i][j].size(); r++)
                    {
                        if ( neighbornode[i][j][r] == vertex[k2] )  //ȷ��k+1������k����ĵڼ����ڵ�
                        {
                            temp1 = Nx[t_l][i][j][r] * (nodev[i][j][0] + nodev[p][q][0]) + Ny[t_l][i][j][r] * (nodev[i][j][1] + nodev[p][q][1]);
                            temp1 = temp1 * L[t_l][i][j][r];

                            temp2 = Nx[t_n][i][j][r] * (nodev[i][j][0] + nodev[p][q][0]) + Ny[t_n][i][j][r] * (nodev[i][j][1] + nodev[p][q][1]);
                            temp2 = temp2 * L[t_n][i][j][r];

                            temp = temp + temp1 + temp2;

                            break;
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
                    i = vertex[k] / (ny + 1);
                    j = vertex[k] % (ny + 1);

                    //����L_{r,r+1}^n * P_{r,r+1/2}^{*,i} * N_{r,r+1}^n
                    {
                        k2 = k + 1;
                        if (k2 > 3)
                        {
                            k2 = k2 - 4;
                        }

                        p = vertex[k2] / (ny + 1);
                        q = vertex[k2] % (ny + 1);

                        for (r=0; r < neighbornode[i][j].size(); r++)
                        {
                            if ( neighbornode[i][j][r] == vertex[k2] )
                            {
                                //ע�������ⷨ������N�ķ����ϵ
                                temp1 = temp1 - L[t_l][i][j][r] * Pstar[i][j][1][r] * Nx[t_l][i][j][r];
                                temp2 = temp2 - L[t_l][i][j][r] * Pstar[i][j][1][r] * Ny[t_l][i][j][r];
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

                        p = vertex[k2] / (ny + 1);
                        q = vertex[k2] % (ny + 1);

                        for (r=0; r<neighbornode[i][j].size(); r++)
                        {
                            if ( neighbornode[i][j][r] == vertex[k2] )
                            {
                                temp1 = temp1 + L[t_l][i][j][r] * Pstar[i][j][0][r] * Nx[t_l][i][j][r];
                                temp2 = temp2 + L[t_l][i][j][r] * Pstar[i][j][0][r] * Ny[t_l][i][j][r];
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
                    i = vertex[k] / (ny + 1);
                    j = vertex[k] % (ny + 1);

                    //�ȼ���L_{r,r+1}^n * P_{r,r+1/2}^{*,i} * V*_r * N_{r,r+1}^n
                    {
                        k2 = k + 1;
                        if (k2 > 3)
                        {
                            k2 = k2 - 4;
                        }

                        p = vertex[k2] / (ny + 1);
                        q = vertex[k2] % (ny + 1);

                        for (r=0; r<neighbornode[i][j].size(); r++)
                        {
                            if ( neighbornode[i][j][r] == vertex[k2] )
                            {
                                temp1 = nodev[i][j][0] * Nx[t_l][i][j][r] + nodev[i][j][1] * Ny[t_l][i][j][r];
                                temp = temp - L[t_l][i][j][r] * Pstar[i][j][1][r] * temp1;  //ע����������ⷨ������N�ķ����ϵ
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

                        p = vertex[k2] / (ny + 1);
                        q = vertex[k2] % (ny + 1);

                        for (r=0; r<neighbornode[i][j].size(); r++)
                        {
                            if ( neighbornode[i][j][r] == vertex[k2] )
                            {
                                temp1 = nodev[i][j][0] * Nx[t_l][i][j][r] + nodev[i][j][1] * Ny[t_l][i][j][r];
                                temp = temp + L[t_l][i][j][r] * Pstar[i][j][0][r] * temp1;
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
