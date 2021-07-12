#include "basic_definitions.h"

double **Allocate_2D(double ** &m, int t1, int t2) {
        m=new double* [t1];
        for (int i=0; i<t1; ++i) {
                m[i]=new double [t2];
                for (int j=0; j<t2; ++j)
                        m[i][j]=0.0;
        }
        return m;
}
