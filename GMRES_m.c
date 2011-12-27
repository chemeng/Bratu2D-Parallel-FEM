//
//  gmres.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <mpi.h>
#include "structs.h"
#include "extern.h"
#include "const.h"
#include "methods_decl.h"


void GMRES_m(double *d,double *r1,struct common4 *sparse,int *indexAA,int rank,int nump,int *inner,MPI_Request *reqs,MPI_Status *stats) 
{
    int i=0,j=0,k=0,l=0,k1=0,k2=0,iter=0,tag1=0,tag2=0,right=0,left=0,ierr;
    double tot_rest=1,vita=0,eps=0,*buffer1,*buffer2,*buffer3,*buffer4,buff=0;
    double e[m+1],y[m],g[m+1];
    double *x,*x0,*r0,*uj,*temp;

    //allocation wste na apothikeutoun sto heap
    x=(double*)calloc(N,sizeof(double));
    x0=(double*)calloc(N,sizeof(double));
    r0=(double*)malloc(N*sizeof(double));
    uj=(double*)malloc(N*sizeof(double));
    temp=(double*)malloc(N*sizeof(double));
    buffer1=(double*)malloc(dim_buf*sizeof(double));
    buffer2=(double*)malloc(dim_buf*sizeof(double));
    buffer3=(double*)malloc(dim_buf*sizeof(double));
    buffer4=(double*)malloc(dim_buf*sizeof(double));
    
    iter=1;
    while (iter<=GMRES_iter) {
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            x0[i]=x[i];
        }
            //midenismos pinakwn-anusmatwn Krylov
        for (i=0; i<N; i++) {
            uj[i]=0;
            r0[i]=0;
            temp[i]=0;
            for (j=0; j<m; j++) {
                gelim.u_base[i][j]=0;
            }
        }
        for (i=0; i<(m+1); i++) {
            e[i]=0;
            g[i]=0;
            for (j=0; j<m; j++) {
                gelim.Hm[i][j]=0;
                y[j]=0;
            }
        }
        e[0]=1;
            //upologismos r0=b-A*x0 opou b=r1 me matmul se CSR format parallel
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            k1=(*sparse).IA[i];
            k2=(*sparse).IA[i+1]-1;
            for (j=k1; j<=k2; j++) {
                r0[i]+=(*sparse).AA[j]*x0[(*sparse).JA[j]];
            }
        }
            //kathe process exei to kommati tou r0 pou upologise+olo t sunoriako
        if (nump>1) {
            for (i=0; i<dim_buf; i++) {
                buffer1[i]=0;
                buffer2[i]=0;
                buffer3[i]=0;
                buffer4[i]=0;
            }
            tag1=11;
            tag2=21;
            left=rank-1;
            right=rank+1;
                //to kathe process exei to kommati tou r0 pou upologise
                //+olokliro sta sunoriaka
            k=0;
            l=0;
            for (i=indexAA[0]; i<=indexAA[1]; i++) {
                if (inner[i]==1) {
                    buffer2[k]=r0[i];
                    k++;
                }
                else if(inner[i]==2){
                    buffer1[l]=r0[i];
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
                    r0[i]+=buffer4[k];
                    k++;
                }
                else if (inner[i]==2){
                    r0[i]+=buffer3[l];
                    l++;
                }
            }
        }
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            r0[i]=r1[i]-r0[i];
        }
            //end CSR multi
        vita=0;
            //upologismos normas tou r0 sto vita->prosoxi sta sunoriaka
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            if ((inner[i]==0) || (inner[i]==1)) {
                vita+=pow(r0[i],2);
            }                  
        }
        buff=0;
        MPI_Allreduce(&vita,&buff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        vita=sqrt(buff);
            //upologismos arxikou dianusmatos uj k apothikeusi ston u+base
        for (i=indexAA[0];i<=indexAA[1]; i++) {
            uj[i]=r0[i]/vita;
        }
            //apothikeusi tou u1 ston pinaka vasewn tou Krylov
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            gelim.u_base[i][0]=uj[i];
        }
            //dianusma g=b*e1
        for (i=0; i<(m+1); i++) {
            g[i]=vita*e[i];
        }
            //kalesma synartisewn
        Krylov(indexAA,inner,sparse,uj,rank,nump,reqs,stats);
        Least_sq(g);
        Back_sub(y,g);
            //x=x0+MATMUL(u_base(N,m),y(m,1)))
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            for (k=0; k<m; k++) {
                temp[i]+=gelim.u_base[i][k]*y[k];
            }
        }
        for (i=indexAA[0]; i<=indexAA[1]; i++) {
            x[i]=x0[i]+temp[i];
        }
        iter++;
    }
    for (i=indexAA[0]; i<=indexAA[1]; i++) {
        d[i]=x[i];
    }
    free(x);
    free(x0);
    free(r0);
    free(uj);
    free(temp);
}
                                                          
                                                          