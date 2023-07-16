#pragma once
#include <string.h>
#include <direct.h>
#include <iostream>
using namespace std;
#define e 2.71828182845904523536
struct para
{
    string rote;
    int row,col;
    int BEnum;//给定的间隔网格数
    double Rside;//右侧极限位置
    double Lside;//左侧极限位置
    double FlyBe;//机翼开始位置
    double H1b,H1e;//水平方向放缩参数
    double H2b,H2e;
    double HenStart;
    double center;
    double H3b,H3e;
    double MaxIteration;
    int Tole;
    para() {
        row = 0;
        col = 0;
        BEnum = 0;
        Rside = 0;
        Lside = 0;
        FlyBe = 0;
        center = 0;
        Tole = 5;
        HenStart = 0;
        rote = "";
        
        H1b = 0;
        H1e = 0;
        H2b = 0;
        H2e = 0;
        H3b = 0;
        H3e = 0;
        MaxIteration = 0;
    }
};

