//
//  main.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
    //Bratu 2D Parallel C implementation
    //NTUA, School of Chemical Engineering
    
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <mpi.h>
/*#include "abfind.c"
#include "Back_sub.c"
#include "CSR.c"
#include "GMRES_m.c"
#include "Krylov.c"
#include "Least_sq.c"
#include "nodnumb.c"
#include "tsfun.c"
#include "xycoord.c"
#include "xydiscr.c"*/
#include "extern.h"
#include "structs.h"
#include "const.h"
#include "methods_decl.h"
    
    //global variable N->unknowns number definition,dim_buf->mpi_buffer size
    int N=0,dim_buf=0,dim=0,Nz=0;
    //global
    struct common3 gelim;

int main (int argc,char* argv[])
{
    int indexAA[2]={0,0},*inner;
    //MPI variables declaration
    MPI_Status stats[4];
    MPI_Request reqs[4];
    double *buffer1,*buffer2,*buffer3,*buffer4;
    int nell=0,i=0,j=0,p=0,q=0,iterat=0,k=0,l=0,init_el=0,fin_el=0;
//    int ncod[nnmax];
    int *ncod;
    int nump=0,rank=0,tag1=0,tag2=0,right=0,left=0;
    double wtime_i=0,wtime_f=0,wtime=0,time=0,tem=0,tem1=0;
    //float xpt[nnmax],ypt[nnmax];
    float *xpt,*ypt;
    int temp_div=0,temp_mod=0;
//    double r1[nnmax],u[nnmax],d[nnmax]
    double *r1,*u,*d;
    struct common1 mesh1;
    struct common2 elem1;
    struct common4 sparse;
    struct common5 param;
    
        //MPI_world build
    MPI_Init( &argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nump);
    wtime_i=MPI_Wtime();
    //dynamic allocation gia na min exw stack-overflow->heap
    //calloc kanei kai initialization me 0
    xpt=(float*)calloc(nnmax,sizeof(double));
    ypt=(float*)calloc(nnmax,sizeof(double));
    ncod=(int*)calloc(nnmax,sizeof(int));
    r1=(double*)calloc(nnmax,sizeof(double));
    u=(double*)calloc(nnmax,sizeof(double));
    d=(double*)calloc(nnmax,sizeof(double));
        //Parametroi provlimatos
    param.alfa=0;
    param.lamda=2;
    xydiscr(&mesh1);
    nodnumb(&elem1);
    xycoord(&mesh1,&elem1,xpt,ypt);
    if (nex<nump) {
        printf("\nCondition nex>nump not met. Parallel solving not possible. Program will now exit.")   ;
        MPI_Finalize();
    }
    N=elem1.np;
    //orismos nzeros-sparsity
    sparse.nzeros=(long int)(0.01*N*(N/nump));
    dim_buf=elem1.nny+2;
    dim=N+1;
       //allocation of dynamic arrays
    sparse.AA=(double*)malloc(sparse.nzeros*sizeof(double));
    sparse.JA=(int*)malloc(sparse.nzeros*sizeof(int));
    sparse.IA=(int*)malloc(dim*sizeof(int));
    inner=(int*)malloc(N*sizeof(int));
    buffer1=(double*)malloc(dim_buf*sizeof(double));
    buffer2=(double*)malloc(dim_buf*sizeof(double));
    buffer3=(double*)malloc(dim_buf*sizeof(double));
    buffer4=(double*)malloc(dim_buf*sizeof(double));
    for (i=0; i<N; i++) {
        inner[i]=(-1);
    }
        //Prepare for Dirichlet boundary conditions
        //markarisma aristera
    for (i=0;i<elem1.nny;i++){
         ncod[i]=1;
    }
        //markarisma katw
    for (i=0;i<(N-elem1.nny+1);(i+=elem1.nny)){
        ncod[i]=1;
    }
        //markarisma panw
    for (i=(elem1.nny-1);i<N;i+=elem1.nny){
        ncod[i]=1;
    }
       //markarisma dexia
    for (i=(N-elem1.nny);i<N;i++){
        ncod[i]=1;
    }
    //orismos kommatiwn elements gia kathe process
    temp_mod=(elem1.ne)%nump;
    temp_div=(elem1.ne)/nump;
    if (temp_mod==0){
        init_el=rank*temp_div;
        fin_el=(rank+1)*temp_div-1;
    }
    else if (temp_mod!=0){
        if (rank!=(nump-1)){
            init_el=(rank*temp_div+rank);
            fin_el=((rank+1)*(temp_div+1)-1);
        }
        else if (rank==(nump-1)){
            init_el=(rank*temp_div+rank);
            fin_el=(elem1.ne-1);
        }
    }
    //indexAA periexei prwti-teleutaia mi mideniki grammi sto kathe process
    indexAA[0]=elem1.nop[init_el][0];
    indexAA[1]=elem1.nop[fin_el][8];
        //ksekinima NEWTON
        iterat=1;
    while (iterat<=Newton_iter){
            //arxikopoiisi prwti fora
        if (iterat==1) {
            for (i=0;i<sparse.nzeros;i++){
                sparse.AA[i]=0;
                sparse.JA[i]=0;
            }
        }
        else {
            for (i=sparse.IA[indexAA[0]];i<=sparse.IA[indexAA[1]+1];i++){
                sparse.AA[i]=0;
                sparse.JA[i]=0;
            }
        }
        for (i=0;i<dim;i++){
            sparse.IA[i]=i;
        }
        for (i=0;i<N;i++){
            r1[i]=0;
            d[i]=0;
        }
        Nz=0;
            //kathe process sarwnei ena kommati tou mesh->sugekrimena elements
        for (nell=init_el;nell<=fin_el;nell++){
            abfind(nell,r1,u,&elem1,xpt,ypt,&sparse,&param,ncod,indexAA);
        }
            //indexAA periexei tin prwti kai teleutaia mi mideniki grammi gia ton AAlocal
        if(iterat==1){
            /* markarisma sunoriakwn nodes me
             Inner nodes=0
             Left boundary=1
             Right boundary =2
             oi times autes mpainoun sto inner */
            if (rank==0) {
                for (i=0;i<elem1.nop[fin_el][5];i++)
                    inner[i]=0;
            }
            if(rank==(nump-1)) {
                for (i=elem1.nop[init_el][3];i<N;i++)
                    inner[i]=0;
            }
            else {
                for (i=elem1.nop[init_el][3];i<=elem1.nop[fin_el][5];i++){
                    inner[i]=0;
                }
            }
            temp_mod=init_el%ney;
                //markarisma aristera nodes me 1
            if (rank!=0) {                
                for (i=0; i<(2*(ney-temp_mod)+1); i++) {
                    inner[elem1.nop[init_el][0]+i]=1;
                }
                if (temp_mod!=0) {
                    inner[elem1.nop[init_el][3]]=1;
                    for (i=0;i<(2*temp_mod+1);i++) {
                        inner[elem1.nop[init_el+ney-temp_mod][0]+i]=1;
                    }
                }
            }
            temp_mod=(fin_el+1)%ney;
            if (rank!=(nump-1)) {
                    //markarisma sunoriakwn dexia nodes me 2
                if (temp_mod==0) {
                    for (i=0; i<(2*ney+1); i++) {
                        inner[elem1.nop[fin_el][8]-i]=2; 
                    }
                }
                if (temp_mod!=0) {
                    for (i=0; i<(2*temp_mod+1); i++) {
                        inner[elem1.nop[fin_el][8]-i]=2;
                    }
                    inner[elem1.nop[fin_el][5]]=2;
                }
                for (i=0; i<(2*(ney-temp_mod)+1); i++) {
                    inner[elem1.nop[fin_el-temp_mod][8]-i]=2;
                }
            }
        }
        if (nump>1) {
            for (i=0; i<dim_buf; i++) {
                buffer1[i]=0;
                buffer2[i]=0;
                buffer3[i]=0;
                buffer4[i]=0;
            }
            tag1=10;
            tag2=10;
            left=rank-1;
            right=rank+1;
                //to kathe process exei to kommati tou r1 pou upologise
                //+olokliro sta sunoriaka
            k=0;
            l=0;
            for (i=indexAA[0]; i<=indexAA[1]; i++) {
                if (inner[i]==1) {
                    buffer2[k]=r1[i];
                    k++;
                }
                else if(inner[i]==2){
                    buffer1[l]=r1[i];
                    l++;
                }
            }
            if (rank==0) {
                MPI_Irecv(buffer3,dim_buf,MPI_DOUBLE,right,tag1,MPI_COMM_WORLD,&reqs[0]);
                MPI_Issend(buffer1,dim_buf,MPI_DOUBLE,right,tag2,MPI_COMM_WORLD,&reqs[1]);
            }
            else if(rank==nump-1){
                MPI_Irecv(buffer4,dim_buf,MPI_DOUBLE,left,tag2,MPI_COMM_WORLD,&reqs[0]);
                MPI_Issend(buffer2,dim_buf,MPI_DOUBLE,left,tag1,MPI_COMM_WORLD,&reqs[1]);
            }
            else {
                MPI_Irecv(buffer4,dim_buf,MPI_DOUBLE,left,tag2,MPI_COMM_WORLD,&reqs[0]);
                MPI_Issend(buffer2,dim_buf,MPI_DOUBLE,left,tag1,MPI_COMM_WORLD,&reqs[1]);
                MPI_Irecv(buffer3,dim_buf,MPI_DOUBLE,right,tag1,MPI_COMM_WORLD,&reqs[2]);
                MPI_Issend(buffer1,dim_buf,MPI_DOUBLE,right,tag2,MPI_COMM_WORLD,&reqs[3]);
            }	
            for (i=0; i<2; i++) {
                MPI_Wait(&reqs[i],&stats[i]);
            }
            if ((rank>0) && (rank<(nump-1))) {
                for (i=2; i<4; i++) {
                    MPI_Wait(&reqs[i],&stats[i]);
                }
            }
            k=0;
            l=0;
            for (i=indexAA[0]; i<=indexAA[1]; i++) {
                if (inner[i]==1) {
                    r1[i]+=buffer4[k];
                    k++;
                }
                else if (inner[i]==2){
                    r1[i]+=buffer3[l];
                    l++;
                }
            }
        }
            //Impose essential boundary conditions gia ton pinaka r1
            //gia ton sk einai mesa stin CSR
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            if (ncod[i]==1) {
                r1[i]=0;
                r1[i]=-(u[i]-0);
            }
        }
        GMRES_m(d,r1,&sparse,indexAA,rank,nump,inner,reqs,stats);
            //kathe process exei to kommati tis lusis tou+olo to sunoriako
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            u[i]+=d[i];
        }
        iterat++;
    }
        //end while
    wtime_f=MPI_Wtime();
    wtime=wtime_f-wtime_i;
    
    if (rank==0) {
        printf("##################################################\n");
        printf("##################################################\n");
        printf("\t\t Parallel Solving of 2-D non linear Bratu Problem\n");
        printf("\n Solver: GMRES Iterative\tMatrix format:CSR\tParallel Library:MPI\n");
        printf("\nWriting process:%d\tout of %d\n",rank,nump);
        printf("\n nex=%d ney=%d ne=%d N=%d\n",nex,ney,elem1.ne,N);
        printf("\n nzeros=%ld, Nz for process zero=%d\n",sparse.nzeros,Nz);
        printf("\n Problem Solved, Newton Converged\n");
    }
        MPI_Allreduce(&wtime,&time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if (rank==0) {
        printf("\nMaximum Time=%lfs\n",time);
        printf("\nNUMP=%d\n",nump);
    }
    if (u[(N+1)/2]!=0) {
        printf("\n#############  SOLUTION  #############");
        printf("\n\n\n u=%.10lf\t lamda=%lf \trank=%d \n\n\n",u[(N+1)/2],param.lamda,rank);
    }
    MPI_Finalize();
return 0;
}

