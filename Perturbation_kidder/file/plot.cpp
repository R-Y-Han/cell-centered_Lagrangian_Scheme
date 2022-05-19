#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

int main()
{

    vector<double > t_tau;
    vector<double > npr;
    vector<double > pr;

    fstream np, p;
    const char* rnp = "cell-centered_Lagrangian_Scheme\\Perturbation_kidder\\output\\r_tn.csv";
    const char* rp = "cell-centered_Lagrangian_Scheme\\Perturbation_kidder\\output\\r_t.csv";

    np.open(rnp, ios :: in);
    p.open(rp, ios :: in);

    string line;
    while (getline(np, line))//getline(inFile, line)表示按行读取CSV文件中的数据
    {
        istringstream sin(line);
        string info;
        vector<string> nnpp;
        while (getline(sin, info, ',')) {
            nnpp.push_back(info);
        }
        
        string np_str = nnpp[1];

        double r;
        stringstream sr;
        sr << np_str;
        sr >> r;
        npr.push_back(r);
    }
    np.close();
    
    while (getline(p, line))//getline(inFile, line)表示按行读取CSV文件中的数据
    {
        istringstream sin(line);
        string info;
        vector<string> pp;
        while (getline(sin, info, ',')) {
            pp.push_back(info);
        }
        
        string tau = pp[0];
        string p_str = pp[1];

        double r,t;
        stringstream sr,st;
        sr << p_str;
        st << tau;
        sr >> r;
        st >> t;

        t_tau.push_back(t);
        pr.push_back(r);
    }
    p.close();

    const char* a = "cell-centered_Lagrangian_Scheme\\Perturbation_kidder\\output\\a.plt";
    fstream a_a0;
    remove(a);
    
    a_a0.open(a,ios:: out | ios :: app);
    for (int h=0; h<t_tau.size(); h++)
    {
        double at;
        at = pr[h] - npr[h];
        at = at / 1e-6;
        a_a0<<t_tau[h]<<"\t"<<at<<endl;
    }
    a_a0.close();

    system("pause");
}