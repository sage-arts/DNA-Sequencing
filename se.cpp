#include<bits/stdc++.h>
using namespace std;

const int N=4380; 
const int M=4469;
const int n=530;

int X[N][M];
int y[N];

int x0[n][M];
int x1[n][M];
int x3[n][M];
int x4[n][M];
int x6[n][M];

int x0t[M][n];
int x1t[M][n];
int x3t[M][n];
int x4t[M][n];
int x6t[M][n];

double cov0[M][M];
double cov1[M][M];
double cov3[M][M];
double cov4[M][M];
double cov6[M][M];

double icov0[M][M];
double icov1[M][M];
double icov3[M][M];
double icov4[M][M];
double icov6[M][M];

double u0[M];
double u1[M];
double u3[M];
double u4[M];
double u6[M];

int vec[M];
double vec1[M];

void mean(int x[][M],double u[])
{
     int i,j;
     for(j=0;j<M;j++)
     {
          for(i=0;i<n;i++)
          {
               u[j]+=x[i][j];
          }
     }
     for(j=0;j<M;j++)
     {
          u[j]/=n;
     }
}

void shift(int x[n][M],double u[])
{
     int i,j;
     for(i=0;i<n;i++)
     {
          for(j=0;j<M;j++)
          {
               x[i][j]-=u[j];
          }
     }
}

void transpose(int x[n][M],int xt[M][n])
{
     int i,j;
     for(i=0;i<n;i++)
     {
          for(j=0;j<M;j++)
          {
               xt[j][i]=x[i][j];
          }
     }   
}

void covar(int x[n][M],int xt[M][n],double cov[M][M])
{
     int i,j,k;
     for(i=0;i<M;i++)
     {
          for(j=0;j<n;j++)
          {
               for(k=0;k<n;k++)
               {
                    cov[i][j]+=xt[i][k]*x[k][j];
               }
          }
     }
     for(i=0;i<M;i++)
     {
          for(j=0;j<n;j++)
          {
               cov[i][j]/=n-1;
          }
     }
}

void transpose1(double cov[M][M],double icov[M][M])
{
     int i,j;
     for(i=0;i<M;i++)
     {
          for(j=0;j<M;j++)
          {
               icov[j][i]=cov[i][j];
          }
     }
}

void dot(int vec[M],double icov[M][M],double vec1[M])
{
     int i,j;
     for(i=0;i<M;i++)
     {
          for(j=0;j<M;j++)
          {
               vec1[i]+=icov[i][j]*vec[j];
          }
     }
}

double dot1(int vec[M],double vec1[M])
{
     double ans=0;
     int i;
     for(i=0;i<M;i++)
     {
          ans+=vec1[i]*vec[i];
     }
     return ans;
}

void clearvec1()
{
     for(int i=0;i<M;i++)
     {
          vec1[i]=0;
     }
}

int main()
{
     string fname="data.csv";
     string line,word;
     fstream file(fname,ios::in);
     if(file.is_open())
     {
          int i=-1;
          while(getline(file,line))
          {
               stringstream str(line);
               int j=-1;
               while(getline(str,word,','))
               {
                    stringstream container(word);
                    if(j==M)
                    {
                         container>>y[i];
                    }
                    else if(i>-1&&j>-1)
                    {
                         container>>X[i][j];
                    }
                    j++;
               }
               i++;
          }
     }

     int k0=0,k1=0,k3=0,k4=0,k6=0;
     for(int i=0;i<N;i++)
     {
          if(y[i]==0)
          {
               for(int j=0;j<M;j++)
               {
                    x0[k0][j]=X[i][j];
               }   
               k0++; 
          }
          else if(y[i]==1)
          {
               for(int j=0;j<M;j++)
               {
                    x1[k1][j]=X[i][j];
               }   
               k1++; 
          }
          else if(y[i]==3)
          {
               for(int j=0;j<M;j++)
               {
                    x3[k3][j]=X[i][j];
               }   
               k3++; 
          }
          else if(y[i]==4)
          {
               for(int j=0;j<M;j++)
               {
                    x4[k4][j]=X[i][j];
               }   
               k4++; 
          }
          else if(y[i]==6)
          {
               for(int j=0;j<M;j++)
               {
                    x6[k6][j]=X[i][j];
               }   
               k6++; 
          }
     }

     mean(x0,u0);
     mean(x1,u1);
     mean(x3,u3);
     mean(x4,u4);
     mean(x6,u6);

     shift(x0,u0);
     shift(x1,u1);
     shift(x3,u3);
     shift(x4,u4);
     shift(x6,u6);

     transpose(x0,x0t);
     transpose(x1,x1t);
     transpose(x3,x3t);
     transpose(x4,x4t);
     transpose(x6,x6t);

     covar(x0,x0t,cov0);
     covar(x1,x1t,cov1);
     covar(x3,x3t,cov3);
     covar(x4,x4t,cov4);
     covar(x6,x6t,cov6);

     transpose1(cov0,icov0);
     transpose1(cov1,icov1);
     transpose1(cov3,icov3);
     transpose1(cov4,icov4);
     transpose1(cov6,icov6);

     double Max=DBL_MIN;
     int label=-1;
     for(int i=0;i<M;i++)
     {
          vec[i]=x6[0][i];
     }
     clearvec1();
     dot(vec,icov0,vec1);
     if(dot1(vec,vec1)>Max) label=0;

     clearvec1();
     dot(vec,icov1,vec1);
     if(dot1(vec,vec1)>Max) label=1;

     clearvec1();
     dot(vec,icov3,vec1);
     if(dot1(vec,vec1)>Max) label=3;

     clearvec1();
     dot(vec,icov4,vec1);
     if(dot1(vec,vec1)>Max) label=4;

     clearvec1();
     dot(vec,icov6,vec1);
     if(dot1(vec,vec1)>Max) label=6;
     cout<<"Class label of the test sequence is: "<<label<<endl;

     return 0;
}