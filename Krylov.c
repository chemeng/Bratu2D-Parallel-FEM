//
//  Krylov.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/20/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <mpi.h>
#include "structs.h"
#include "const.h"
#include "extern.h"
#include "methods_decl.h"

void Krylov(int *indexAA,int *inner,struct common4 *sparse,double *uj,int rank,int nump,MPI_Request *reqs,MPI_Status *stats)
{
    int  i=0,j=0,k=0,k1=0,k2=0,p=0,tag1=0,tag2=0,right=0,left=0,count=0,ierr,l=0;   
    double *buffer1,*buffer2,*buffer3,*buffer4,buff=0,buff1=0;
    double *w;
    
    buffer1=(double*)malloc(dim_buf*sizeof(double));
    buffer2=(double*)malloc(dim_buf*sizeof(double));
    buffer3=(double*)malloc(dim_buf*sizeof(double));
    buffer4=(double*)malloc(dim_buf*sizeof(double));
    
    w=(double*)malloc(N*sizeof(double));
    for (j=0;j<m;j++){
            //j einai to count g mas
            //vazei to kommati tou uj pou antistoixei sto process
        for(k=indexAA[0];k<=indexAA[1];k++){
            if(inner[k]!=(-1)){
                uj[k]=gelim.u_base[k][j];
            }
        }
            //matmul me CSR format w=MATMUL(A,uj)
        for(i=0;i<N;i++){
            w[i]=0;
        }
        for(i=indexAA[0];i<=indexAA[1];i++){
            if(inner[i]!=(-1)){
                k1=(*sparse).IA[i];
                k2=(*sparse).IA[i+1];
                for (k=k1;k<k2;k++){
                    w[i]+=((*sparse).AA[k])*(uj[(*sparse).JA[k]]);
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
            tag1=31;
            tag2=41;
            left=rank-1;
            right=rank+1;
                //to kathe process exei to kommati tou w pou upologise
                //+olokliro sta sunoriaka
            k=0;
            l=0;
            for (i=indexAA[0]; i<=indexAA[1]; i++) {
                if (inner[i]==1) {
                    buffer2[k]=w[i];
                    k++;
                }
                else if(inner[i]==2){
                    buffer1[l]=w[i];
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
                    w[i]+=buffer4[k];
                    k++;
                }
                else if (inner[i]==2){
                    w[i]+=buffer3[l];
                    l++;
                }
            }
        }
            //end CSR multi
        //MEXRI EDW SWSTES KAI OI J epanalipseis
        for (i=0;i<=j;i++){
            for(k=indexAA[0];k<=indexAA[1];k++){
                if (inner[k]!=(-1)){
                    uj[k]=gelim.u_base[k][i];
                }
            }
                //DOT_PRODUCT(w,uj)
            gelim.Hm[i][j]=0;
            for (k=indexAA[0];k<=indexAA[1];k++){
                if((inner[k]==0) || (inner[k]==1)){
                    gelim.Hm[i][j]+=w[k]*uj[k];
                }
            }
            buff=0;
            buff1=gelim.Hm[i][j];
            MPI_Allreduce(&buff1,&buff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            gelim.Hm[i][j]=buff;
                //end DOT product
            for(k=indexAA[0];k<=indexAA[1];k++){
                if (inner[k]!=-1){
                    w[k]-=(gelim.Hm[i][j])*(uj[k]);
                }
            }
        }
        for (k=indexAA[0];k<=indexAA[1];k++){
            if((inner[k]==0) || (inner[k]==1)){
                gelim.Hm[j+1][j]+=pow(w[k],2);
            }
        }
        buff=0;
        buff1=gelim.Hm[j+1][j];
        MPI_Allreduce(&buff1,&buff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        gelim.Hm[j+1][j]=sqrt(buff);
        if(j<(m-1)){
            for (k=indexAA[0];k<=indexAA[1];k++){
                if(inner[k]!=-1){
                    gelim.u_base[k][j+1]=w[k]/(gelim.Hm[j+1][j]);
                }
            }
        } 
    }
}   

