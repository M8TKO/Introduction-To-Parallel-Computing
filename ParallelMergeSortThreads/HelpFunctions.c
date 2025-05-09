#include "HelpFunctions.h"
#include <stdio.h>
int binary_search(double *a, int p, int r, double x){
    int low = p, high = r+1, mid, q = p;
    while( low < high ){

        mid = ( low + high ) / 2;
        if( x <= a[mid] )
            high = mid;
        else    
            low = mid + 1;
        q = low;
    } 
    return q;
}

void sequential_merge(double *a, int p1,int r1, int p2, int r2, double *b, int p3){
    int n1 = r1 - p1 + 1,  
        n2 = r2 - p2 + 1,
        i = 0, j = 0, k = p3;
    
    while( i < n1 && j < n2 ){
        
        if( a[p1 + i] <= a[p2 + j] )
            b[k] = a[p1 + (i++)];
        else    
            b[k] = a[p2 + (j++)];
        k++;
    }

    while ( i < n1 ){
        b[k++] = a[p1 + (i++)];
    }
    while ( j < n2 ){
        b[k++] = a[p2 + (j++)];
    }
}