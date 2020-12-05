#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    //VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);
    MatrixXd P_xyz = MatrixXd::Zero(m*p_num1d , 3) ;
    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    // 映射矩阵A为单个维度矩阵，三个维度相同
    MatrixXd A = MatrixXd::Zero(m * p_num1d, m * p_num1d);
    VectorXd time_pow = VectorXd::Zero(p_num1d);
    VectorXd fac_table = VectorXd::Zero(p_num1d);//Factorial查找表
    for (int j=0;j < p_num1d;j++)
        {
            fac_table(j) = Factorial(j);
        } 
    //cout<<"Factorial查找表: \n"<<fac_table.transpose()<<endl;
    MatrixXd row_A = MatrixXd::Zero(p_num1d , p_num1d);
    //cout<<"开始规划：\n";
    for (int i =0;i < m;i++)
    {
        double time = Time(i);
        //cout<<"time : "<<time<<endl;
        //row_A上半部分
        for (int j=0;j < d_order;j++)
        {
            row_A(j , p_num1d - j -1) = fac_table(j);
        } 
        //row_A下半部分
        //cout<<"row_A : \n"<<row_A<<endl;
        for (int i =0; i < p_num1d;i++)//计算时间乘方
            time_pow(i) = pow(time, i);
        //cout<<"time_pow: \n"<<time_pow.transpose()<<endl;
        for (int i=0;i<p_num1d;i++)
        {
            for (int j =0;j< d_order;j++)
                if (p_order - i -j >= 0)
                row_A(j + d_order, i) = fac_table(p_order - i)/fac_table(p_order - i -j)
                                                 * time_pow(p_order - i -j);
        }
        //cout<<"row_A: \n"<<row_A<<endl;
        A.block(i* p_num1d, i*p_num1d, p_num1d,p_num1d) = row_A;
    }
    //cout<<"映射矩阵A： \n"<<A<<endl;
    //计算Q矩阵，单个维度，三个维度相同
    MatrixXd Q = MatrixXd::Zero(m * p_num1d, m * p_num1d);
    MatrixXd _Q = MatrixXd::Zero(p_num1d, p_num1d);//临时存储
    for (int i=0;i<m;i++)
    {
        for (int j=0;j<p_num1d;j++)
        {
            for (int k=0;k<p_num1d;k++)
            {
                if(j >= d_order && k >= d_order)
                {
                    _Q(p_num1d -1 -j, p_num1d -1 -k) = (fac_table(j)*fac_table(k)*time_pow(j + k - p_order))/
                                (fac_table(j- d_order) * fac_table(k - d_order) * (j + k - p_order));
                }
            }
        }
        //cout<<"_Q: \n"<<_Q<<endl;
        Q.block(i* p_num1d, i*p_num1d,p_num1d,p_num1d) = _Q;
    }
    //cout<<"系数矩阵Q为： \n"<<Q<<endl;
    //计算C矩阵，三个维度独立，仅计算一个维度
    MatrixXd C;//C矩阵
    MatrixXd C_T = MatrixXd::Zero(2 * m * d_order, (m + 1) * d_order);//C转置
    int df_num = 2 * d_order + (m - 1) ;//起点终点位置、速度、加速度约束以及中间节点位置约束
    int dp_num = (m + 1 ) * d_order - df_num;//自由变量数目
    C_T.block(0,0,d_order,d_order) = MatrixXd::Identity(d_order, d_order);//初始约束
    C_T.block((2 * m - 1) * d_order, df_num - d_order, d_order, d_order) = 
                                MatrixXd::Identity(d_order, d_order);//终点约束
    int df_index = d_order;
    int dp_index = df_num; 
    //cout<<"m: "<<m<<"d_orde : "<<d_order<<endl;
    //中间点                           
    for (int i=1; i < 2 * m -1  ; i += 2)
    {     
        C_T(i* d_order , df_index) = 1;//位置约束
        C_T((i + 1 ) * d_order , df_index) = 1;//位置约束
        for (int j=1;j<d_order;j++)
        {
            C_T(i * d_order + j , dp_index) = 1;
            C_T((i+1) * d_order + j , dp_index) = 1;
            dp_index += 1;
        }
        df_index += 1;
    }
    
    C = C_T.transpose();
    //cout<<"选择矩阵C转置:\n "<<C_T<<endl;
    /*   Produce the dereivatives in X, Y and Z axis directly.  */
    //三轴独立，仅求解单个维度

    MatrixXd A_inv = A.inverse();
    MatrixXd A_inv_T = A_inv.transpose();
    MatrixXd R = C * A_inv_T * Q * A_inv * C_T;

    MatrixXd R_fp = R.block(0 , df_num, df_num, dp_num);
    MatrixXd R_pp = R.block(df_num , df_num, dp_num, dp_num);
    
    VectorXd df = VectorXd::Zero(df_num);//约束条件
    VectorXd dp = VectorXd::Zero(dp_num);//优化变量
    VectorXd dfp = VectorXd::Zero(dp_num + df_num);//全部变量，用于计算多项式系数矩阵
    VectorXd state_temp = VectorXd::Zero(d_order);//临时存储变量
    for (int i=0;i<3;i++)
    {
        state_temp(0) = Path.row(0)(i);//起始点
        state_temp(1) = Vel(0,i);
        state_temp(2) = Acc(0,i);
        df.head(d_order) = state_temp;
        df.segment(d_order, m - 1) = Path.col(i).segment(1, m-1);
        state_temp(0) = Path.row(m)(i);//终止点
        state_temp(1) = Vel(1,i);
        state_temp(2) = Acc(1,i);
        df.tail(d_order) = state_temp;
        //cout<<"df: \n"<<df.transpose()<<endl;
        dp = -R_pp.inverse() * R_fp.transpose() * df;
        dfp<<df,dp;
        P_xyz.col(i) = A_inv * C_T * dfp;
    }
    for (int i=0;i<m;i++)
    {
        PolyCoeff.row(i).segment(0, p_num1d) = P_xyz.col(0).segment(i * p_num1d,p_num1d);
        PolyCoeff.row(i).segment(p_num1d, p_num1d) = P_xyz.col(1).segment(i * p_num1d,p_num1d);
        PolyCoeff.row(i).segment(2 * p_num1d, p_num1d) = P_xyz.col(2).segment(i * p_num1d,p_num1d);
    }
    
    return PolyCoeff;
}
