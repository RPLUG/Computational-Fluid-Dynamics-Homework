#include <iostream>
#include "GeneF.h"
#include "DataS.h"

int main()
{
    string rote = GetCurrentPath();
    para para;
    para.row =100;//下面放缩参数A1B1按照这个计算
    para.col = 500;//下面,col需要为偶数，再计算处理时候会自动处理为奇数，方便计算，也就是199,200能保证捕捉到翼型
    para.BEnum = 20;
    para.rote = rote;
    para.Rside = 4;
    para.Lside = 0;//左侧极限位置
    para.HenStart = 1;//用来确定初始网格上下水平的起始位置
    para.Tole = 5;
    para.FlyBe = 1;//翼型头部起始位置
    para.H1b = 0.01;
    para.H1e = 0.95;//同时确定尾部网格的长度
    para.H2b = 0.0001;
    para.H2e = 3;
    para.H3e = 4;
    para.H3b = 0.0001;
    para.MaxIteration = 1000;
    para.center = 0;
    MeshGenerate* Solution = new MeshGenerate(para);
    return 0;
}
