#include <iostream>
#include <cmath>
#include "initial.h"

using namespace std;

//仅针对本例中的条件
int main()
{
    double P, al, ald, ar, ard, f, fd, temp, temp2;
    P = 0.8;
    temp = 0.3;
    while ( abs(P - temp)>1e-9)
    {
        temp = P;
        
        al = sqrt(7.0/5.0) * (1-P) / (7 * (1- pow(P,1.0/7.0)));
        temp2 = 8 * (P - 0.1) + 1;
        ar = sqrt(3) * sqrt(temp2)/12.0;

        f = (P-0.1)/ar + (P-1)/al;
        ald = sqrt(7.0/5.0) * (-7 * (1-pow(P,1.0/7.0)) + (1- P ) * pow(P,-6.0/7.0) );
        ard = 1 / sqrt(3 * (8*P - 0.2));
        fd = (ar - (P-0.1)*ard)/(ar*ar) + (al - (P-1)*ald)/(al*al);

        P = temp - f / fd;
        cout<<P<<endl;

    }
    cout<<"P="<<P<<endl;

    double U;
    temp2 = 8 * (P - 0.1) + 1;
    ar = sqrt(3) * sqrt(temp2)/12.0;
    U = (P-0.1)/ar;
    cout<<"U="<<U<<endl;

    double s;
    s = ar/rho_r + 0;
    cout<<"shock speed="<<s<<endl;

    double rhor;
    rhor = ar / (s-U);
    cout<<"rho*r = "<<rhor<<endl;

    double crare;
    crare = sqrt(5.0/3.0) + 0.2 * (0 - U);
    
    double rhol;
    rhol = 5 * P / (3*crare * crare);
    cout<<"rho*l = "<<rhol<<endl;
    cout<<"rarefaction head speed="<< - sqrt(5.0/3.0)<<endl;
    cout<<"rarefaction tail speed="<<U - crare<<endl;

    system("pause");
}