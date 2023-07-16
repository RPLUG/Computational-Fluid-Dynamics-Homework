#include <iostream>
#include "Solution.h"
#include "BasicData.h"
int main(int argc,char* argv[])
{
    Para para;
    para.DeltX = 0.01;
    para.DeltY = 0.01;
    para.End_Time = 0.1;
    para.Ite_Step = 0.00001;
    para.Xnumber = 100;
    para.Ynumber = 100;
    para.Print = true;
    para.Re = 100;
    para.Inter_ite_time = 10;
    para.Cal_P_Flag = false;
    para.out_route = GetCurrentPath()+"/Result";
    Solution_SCF* NewSol = new Solution_SCF(para);
    return 0;
}
