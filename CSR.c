//
//  csr.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <math.h> 
#include <mpi.h>
#include "structs.h"
#include "extern.h"
#include "const.h"
#include "methods_decl.h"


void CSR(struct common4 *sparse,double value,int i,int j,int *ncod,int *indexAA)
{
    int k=0,l=0,t=0,update=0,temp=0;
    temp=(*sparse).IA[i];
        //Dirichlet boundary condition
    if ((ncod[i]==1) && (i!=j)){
        value=0;
    }
    else if ((ncod[i]==1) && (i==j) && ((*sparse).AA[temp]!=1)) {
        value=1;
    }
    else if ((ncod[i]==1)   &&  (i==j)  && ((*sparse).AA[temp]==1)){
        value=0;
    }    
        //ksekinaei i eisagwgi
    if (value!=0) {
        if ((*sparse).AA[temp]==0) {
            (*sparse).AA[temp]=value;
            (*sparse).JA[temp]=j;
            (Nz)++;
        }
        else if((*sparse).AA[temp]!=0){
            t=(*sparse).IA[i];
            while (t<((*sparse).IA[i+1])) {
                if ((*sparse).JA[t]==j) {
                    update=1;
                    break;
                }
                t++;
            }
        }
        if (update==1) {
            (*sparse).AA[t]+=value;
        }
        else {
            if (((*sparse).AA[temp]!=0) && ((*sparse).JA[temp]>j)) {
                for (k=((*sparse).IA[indexAA[1]+1]); k>=temp; k--) {
                    (*sparse).AA[k+1]=(*sparse).AA[k];
                    (*sparse).JA[k+1]=(*sparse).JA[k];
                }
                (*sparse).AA[temp]=value;
                (*sparse).JA[temp]=j;
                for (k=(i+1); k<(indexAA[1]+2); k++) {
                    ((*sparse).JA[k])++;
                }
                (Nz)++;
            }
            else if ((((*sparse).AA[temp])!=0) && (((*sparse).JA[temp])<j)) {
                k=temp;
                while (((*sparse).JA[k]<j) && (k<((*sparse).IA[i+1]))) {
                    k++;
                }
                for (l=((*sparse).IA[indexAA[1]+1]); l>=k; l--) {
                    (*sparse).AA[l+1]=(*sparse).AA[l];
                    (*sparse).JA[l+1]=(*sparse).JA[l];
                }
                (*sparse).AA[k]=value;
                (*sparse).JA[k]=j;
                for (k=(i+1); k<(indexAA[1]+2); k++) {
                    (*sparse).IA[k]++;    
                }
                (Nz)++;
            }
        }
    }
    /*l=(*sparse).IA[dim-2]+1;
    while ((*sparse).AA[l]!=0) {
        l++;    
    }
    (*sparse).IA[dim-1]=l;
    */
}

