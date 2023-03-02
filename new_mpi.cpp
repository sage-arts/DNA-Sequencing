#include<mpi.h>
#include<bits/stdc++.h>
using namespace std;

#define MASTER 0 
#define FROM_MASTER 1 
#define FROM_WORKER 2 

const int N=4380; 
const int M=4469;
const int n=530;

int X[N][M];
int y[N];

int temp0[n][M];
int temp1[n][M];
int temp3[n][M];
int temp4[n][M];
int temp6[n][M];

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

double pu[M];

int vec[M];
double vec1[M];

double pans, ans;

double Max;
int label;

int numtasks,
taskid, 
numworkers, 
source, 
dest, 
mtype, 
rows, 
averow, extra, offset, 
i, j, k, rc; 
MPI_Status status;
double start, finish;

int main(int argc, char *argv[])
{
     MPI_Init(&argc,&argv);
     MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
     MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
     if (numtasks < 2 ){
          printf("Need at least two MPI tasks. Quitting...\n");
          MPI_Abort(MPI_COMM_WORLD, rc);
          exit(1);
     }
     numworkers = numtasks-1;
     if(taskid == MASTER)
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
                         temp0[k0][j]=X[i][j];
                    }   
                    k0++; 
               }
               else if(y[i]==1)
               {
                    for(int j=0;j<M;j++)
                    {
                         temp1[k1][j]=X[i][j];
                    }   
                    k1++; 
               }
               else if(y[i]==3)
               {
                    for(int j=0;j<M;j++)
                    {
                         temp3[k3][j]=X[i][j];
                    }   
                    k3++; 
               }
               else if(y[i]==4)
               {
                    for(int j=0;j<M;j++)
                    {
                         temp4[k4][j]=X[i][j];
                    }   
                    k4++; 
               }
               else if(y[i]==6)
               {
                    for(int j=0;j<M;j++)
                    {
                         temp6[k6][j]=X[i][j];
                    }   
                    k6++; 
               }
          }
          start = MPI_Wtime();
     }

     // mean(temp0,u0)
     if(taskid==MASTER)
     {
          averow = n/numworkers;
          extra = n%numworkers;
          offset = 0;
          mtype = FROM_MASTER;
          for(dest=1; dest<=numworkers; dest++){
          rows = (dest <= extra) ? averow+1 : averow;
          MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&temp0[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
          offset = offset + rows;
          }
          mtype = FROM_WORKER;
          for(i=1; i<=numworkers; i++){
               source = i;
               MPI_Recv(&pu, M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
               for(j=0;j<M;j++){
                    u0[j]+=pu[j];
               }
          }
          for(j=0;j<M;j++){
               u0[j]/=n;
          }
     }
     if(taskid > MASTER)
     {
          mtype = FROM_MASTER;
          MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&temp0, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          memset(pu,0,sizeof(pu));
          for(i=0;i<rows;i++){
               for(j=0;j<M;j++){
                    pu[j] += temp0[i][j];
               }
          }
          mtype = FROM_WORKER;
          MPI_Send(&pu, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     }

     // shift(temp0,x0,u0);
     if(taskid==MASTER)
     {
          averow = n/numworkers;
          extra = n%numworkers;
          offset = 0;
          mtype = FROM_MASTER;
          for(dest=1; dest<=numworkers; dest++){
          rows = (dest <= extra) ? averow+1 : averow;
          MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&temp0[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
          MPI_Send(&u0, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
          offset = offset + rows;
          }
          mtype = FROM_WORKER;
          for(i=1; i<=numworkers; i++){
               source = i;
               MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&x0[offset][0], rows*M, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
          }
     }
     if(taskid > MASTER)
     {
          mtype = FROM_MASTER;
          MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&temp0, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&u0, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
          for(i=0;i<rows;i++){
               for(j=0;j<M;j++){
                    x0[i][j] = temp0[i][j]-u0[j];
               }
          }
          mtype = FROM_WORKER;
          MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&x0, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     }

     // transpose(x0,x0t);
     if(taskid == MASTER)
     {
          averow = M/numworkers;
          extra = M%numworkers;
          offset = 0;
          mtype = FROM_MASTER;
          for(dest=1; dest<=numworkers; dest++){
          rows = (dest <= extra) ? averow+1 : averow;
          MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&x0, n*M, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          offset = offset + rows;
          }
          mtype = FROM_WORKER;
          for(i=1; i<=numworkers; i++){
               source = i;
               MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&x0t[offset][0], rows*n, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
          }
     }
     if(taskid > MASTER)
     {
          mtype = FROM_MASTER;
          MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&x0, n*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          for(i=0;i<rows;i++){
               for(j=0;j<n;j++){
                    x0t[i][j] = x0[j][i+offset];
               }
          }
          mtype = FROM_WORKER;
          MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&x0t, rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     }

     // covar(x0,x0t,cov0);
     if(taskid == MASTER)
     {
          averow = M/numworkers;
          extra = M%numworkers;
          offset = 0;
          mtype = FROM_MASTER;
          for(dest=1; dest<=numworkers; dest++){
          rows = (dest <= extra) ? averow+1 : averow;
          MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&x0t[offset][0], rows*n, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
          MPI_Send(&x0, n*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
          offset = offset + rows;
          }
          mtype = FROM_WORKER;
          for(i=1; i<=numworkers; i++){
               source = i;
               MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&cov0[offset][0], rows*M, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
          }
     }
     if(taskid > MASTER)
     {
          mtype = FROM_MASTER;
          MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&x0t, rows*n, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&x0, n*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
          for(k=0; k<M; k++){
               for(i=0; i<rows; i++){
                    cov0[i][k] = 0.0;
                    for(j=0; j<n; j++){
                         cov0[i][k] += x0t[i][j] * x0[j][k];
                    }
                    cov0[i][k] /= n-1;
               }
          }
          mtype = FROM_WORKER;
          MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&cov0, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     }

     // // transpose1(cov0,icov0);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov0, M*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&icov0[offset][0], rows*M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&cov0, M*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                icov0[i][j] = cov0[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov0, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // mean(temp1,u1)
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp1[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pu, M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           for(j=0;j<M;j++){
     //                u1[j]+=pu[j];
     //           }
     //      }
     //      for(j=0;j<M;j++){
     //           u1[j]/=n;
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp1, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      memset(pu,0,sizeof(pu));
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                pu[j] += temp1[i][j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pu, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // shift(temp1,x1,u1);
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp1[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&u1, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x1[offset][0], rows*M, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp1, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&u1, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                x1[i][j] = temp1[i][j]-u1[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x1, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose(x1,x1t);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x1, n*M, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x1t[offset][0], rows*n, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x1, n*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                x1t[i][j] = x1[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x1t, rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // covar(x1,x1t,cov1);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x1t[offset][0], rows*n, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&x1, n*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&cov1[offset][0], rows*M, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x1t, rows*n, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x1, n*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(k=0; k<M; k++){
     //           for(i=0; i<rows; i++){
     //                cov1[i][k] = 0.0;
     //                for(j=0; j<n; j++){
     //                     cov1[i][k] += x1t[i][j] * x1[j][k];
     //                }
     //                cov1[i][k] /= n-1;
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov1, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose1(cov1,icov1);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov1, M*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&icov1[offset][0], rows*M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&cov1, M*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                icov1[i][j] = cov1[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov1, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // mean(temp3,u3)
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp3[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pu, M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           for(j=0;j<M;j++){
     //                u3[j]+=pu[j];
     //           }
     //      }
     //      for(j=0;j<M;j++){
     //           u3[j]/=n;
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp3, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      memset(pu,0,sizeof(pu));
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                pu[j] += temp3[i][j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pu, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // shift(temp3,x3,u3);
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp3[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&u3, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x3[offset][0], rows*M, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp3, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&u3, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                x3[i][j] = temp3[i][j]-u3[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x3, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose(x3,x3t);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x3, n*M, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x3t[offset][0], rows*n, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x3, n*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                x3t[i][j] = x3[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x3t, rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // covar(x3,x3t,cov3);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x3t[offset][0], rows*n, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&x3, n*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&cov3[offset][0], rows*M, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x3t, rows*n, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x3, n*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(k=0; k<M; k++){
     //           for(i=0; i<rows; i++){
     //                cov3[i][k] = 0.0;
     //                for(j=0; j<n; j++){
     //                     cov3[i][k] += x3t[i][j] * x3[j][k];
     //                }
     //                cov3[i][k] /= n-1;
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov3, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose1(cov3,icov3);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov3, M*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&icov3[offset][0], rows*M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&cov3, M*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                icov3[i][j] = cov3[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov3, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // mean(temp4,u4)
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp4[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pu, M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           for(j=0;j<M;j++){
     //                u4[j]+=pu[j];
     //           }
     //      }
     //      for(j=0;j<M;j++){
     //           u4[j]/=n;
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp4, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      memset(pu,0,sizeof(pu));
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                pu[j] += temp4[i][j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pu, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // shift(temp4,x4,u4);
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp4[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&u4, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x4[offset][0], rows*M, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp4, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&u4, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                x4[i][j] = temp4[i][j]-u4[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x4, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose(x4,x4t);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x4, n*M, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x4t[offset][0], rows*n, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x4, n*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                x4t[i][j] = x4[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x4t, rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // covar(x4,x4t,cov4);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x4t[offset][0], rows*n, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&x4, n*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&cov4[offset][0], rows*M, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x4t, rows*n, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x4, n*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(k=0; k<M; k++){
     //           for(i=0; i<rows; i++){
     //                cov4[i][k] = 0.0;
     //                for(j=0; j<n; j++){
     //                     cov4[i][k] += x4t[i][j] * x4[j][k];
     //                }
     //                cov4[i][k] /= n-1;
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov4, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose1(cov4,icov4);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov4, M*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&icov4[offset][0], rows*M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&cov4, M*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                icov4[i][j] = cov4[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov4, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // mean(temp6,u6)
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp6[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pu, M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           for(j=0;j<M;j++){
     //                u6[j]+=pu[j];
     //           }
     //      }
     //      for(j=0;j<M;j++){
     //           u6[j]/=n;
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp6, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      memset(pu,0,sizeof(pu));
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                pu[j] += temp6[i][j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pu, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // shift(temp6,x6,u6);
     // if(taskid==MASTER)
     // {
     //      averow = n/numworkers;
     //      extra = n%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&temp6[offset][0], rows*M, MPI_INT, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&u6, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x6[offset][0], rows*M, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&temp6, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&u6, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<M;j++){
     //                x6[i][j] = temp6[i][j]-u6[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x6, rows*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose(x6,x6t);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x6, n*M, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&x6t[offset][0], rows*n, MPI_INT, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x6, n*M, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                x6t[i][j] = x6[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x6t, rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // covar(x6,x6t,cov6);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&x6t[offset][0], rows*n, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&x6, n*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&cov6[offset][0], rows*M, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x6t, rows*n, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&x6, n*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(k=0; k<M; k++){
     //           for(i=0; i<rows; i++){
     //                cov6[i][k] = 0.0;
     //                for(j=0; j<n; j++){
     //                     cov6[i][k] += x6t[i][j] * x6[j][k];
     //                }
     //                cov6[i][k] /= n-1;
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov6, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // transpose1(cov6,icov6);
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&cov6, M*M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&icov6[offset][0], rows*M, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&cov6, M*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0;i<rows;i++){
     //           for(j=0;j<n;j++){
     //                icov6[i][j] = cov6[j][i+offset];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov6, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     if(taskid==MASTER)
     {
          Max=DBL_MIN;
          label=-1;
          for(int i=0;i<M;i++)
          {
               vec[i]=x6[0][i];
          }
     }

     // dot(vec,icov0,vec1)
     if(taskid == MASTER)
     {
          averow = M/numworkers;
          extra = M%numworkers;
          offset = 0;
          mtype = FROM_MASTER;
          for(dest=1; dest<=numworkers; dest++){
          rows = (dest <= extra) ? averow+1 : averow;
          MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&icov0[offset][0], rows*M, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
          MPI_Send(&vec, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
          offset = offset + rows;
          }
          mtype = FROM_WORKER;
          for(i=1; i<=numworkers; i++){
               source = i;
               MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
               MPI_Recv(&vec1[offset], rows, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
          }
     }
     if(taskid > MASTER)
     {
          mtype = FROM_MASTER;
          MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&icov0, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&vec, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
          for(i=0; i<rows; i++){
               vec1[i] = 0.0;
               for(j=0; j<M; j++){
                    vec1[i] += icov0[i][j] * vec[j];
               }
          }
          mtype = FROM_WORKER;
          MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
          MPI_Send(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     }

     // dot1(vec,vec1)
     if(taskid==MASTER)
     {
          ans =0;
          averow = M/numworkers;
          extra = M%numworkers;
          offset = 0;
          mtype = FROM_MASTER;
          for(dest=1; dest<=numworkers; dest++){
          rows = (dest <= extra) ? averow+1 : averow;
          MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
          MPI_Send(&vec[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
          MPI_Send(&vec1[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);      
          offset = offset + rows;
          }
          mtype = FROM_WORKER;
          for(i=1; i<=numworkers; i++){
               source = i;
               MPI_Recv(&pans, 1, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
               ans+=pans;
          }
          if(ans>Max) label=0;
     }
     if(taskid > MASTER)
     {
          mtype = FROM_MASTER;
          MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&vec, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
          MPI_Recv(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
          pans=0;
          for(i=0;i<rows;i++){
               pans += vec[i]*vec1[i];
          }
          mtype = FROM_WORKER;
          MPI_Send(&pans, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     }

     // // dot(vec,icov1,vec1)
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov1[offset][0], rows*M, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&vec1[offset], rows, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&icov1, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0; i<rows; i++){
     //           vec1[i] = 0.0;
     //           for(j=0; j<M; j++){
     //                vec1[i] += icov1[i][j] * vec[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // dot1(vec,vec1)
     // if(taskid==MASTER)
     // {
     //      ans =0;
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec1[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);      
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pans, 1, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           ans+=pans;
     //      }
     //      if(ans>Max) label=1;
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      pans=0;
     //      for(i=0;i<rows;i++){
     //           pans += vec[i]*vec1[i];
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pans, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // dot(vec,icov3,vec1)
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov3[offset][0], rows*M, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&vec1[offset], rows, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&icov3, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0; i<rows; i++){
     //           vec1[i] = 0.0;
     //           for(j=0; j<M; j++){
     //                vec1[i] += icov3[i][j] * vec[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // dot1(vec,vec1)
     // if(taskid==MASTER)
     // {
     //      ans =0;
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec1[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);      
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pans, 1, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           ans+=pans;
     //      }
     //      if(ans>Max) label=3;
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      pans=0;
     //      for(i=0;i<rows;i++){
     //           pans += vec[i]*vec1[i];
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pans, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // dot(vec,icov4,vec1)
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov4[offset][0], rows*M, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&vec1[offset], rows, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&icov4, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0; i<rows; i++){
     //           vec1[i] = 0.0;
     //           for(j=0; j<M; j++){
     //                vec1[i] += icov4[i][j] * vec[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // dot1(vec,vec1)
     // if(taskid==MASTER)
     // {
     //      ans =0;
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec1[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);      
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pans, 1, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           ans+=pans;
     //      }
     //      if(ans>Max) label=4;
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      pans=0;
     //      for(i=0;i<rows;i++){
     //           pans += vec[i]*vec1[i];
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pans, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // dot(vec,icov6,vec1)
     // if(taskid == MASTER)
     // {
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&icov6[offset][0], rows*M, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec, M, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
     //           MPI_Recv(&vec1[offset], rows, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //      }
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&icov6, rows*M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, M, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      for(i=0; i<rows; i++){
     //           vec1[i] = 0.0;
     //           for(j=0; j<M; j++){
     //                vec1[i] += icov6[i][j] * vec[j];
     //           }
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     // // dot1(vec,vec1)
     // if(taskid==MASTER)
     // {
     //      ans =0;
     //      averow = M/numworkers;
     //      extra = M%numworkers;
     //      offset = 0;
     //      mtype = FROM_MASTER;
     //      for(dest=1; dest<=numworkers; dest++){
     //      rows = (dest <= extra) ? averow+1 : averow;
     //      MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
     //      MPI_Send(&vec[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
     //      MPI_Send(&vec1[offset], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);      
     //      offset = offset + rows;
     //      }
     //      mtype = FROM_WORKER;
     //      for(i=1; i<=numworkers; i++){
     //           source = i;
     //           MPI_Recv(&pans, 1, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);
     //           ans+=pans;
     //      }
     //      if(ans>Max) label=6;
     // }
     // if(taskid > MASTER)
     // {
     //      mtype = FROM_MASTER;
     //      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      MPI_Recv(&vec1, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
     //      pans=0;
     //      for(i=0;i<rows;i++){
     //           pans += vec[i]*vec1[i];
     //      }
     //      mtype = FROM_WORKER;
     //      MPI_Send(&pans, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
     // }

     if(taskid==MASTER)
     {
          // cout<<"Class label of the test sequence is: "<<label<<endl;
          finish = MPI_Wtime();
          cout<<"Time= "<<finish - start<<endl;
     }
     MPI_Finalize();
     return 0;
}