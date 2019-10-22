package utils;

import java.util.ArrayList;

public class SquareMatrix extends Matrix {

    public SquareMatrix(int m) {
        super(m,m);
    }

    public SquareMatrix(SquareMatrix that) {
        super(that);
    }

    public SquareMatrix(Matrix that) {
        super(that);
    }

    public SquareMatrix(double[][] data) {
        super(data);

    }

    public SquareMatrix(int m, boolean identity) {
        super(m,m);
        if (identity) {
            for(int i = 0; i < m ; i++) {
                A[i][i] = 1;
            }
        }
    }

    /**
     * Determinant
     * @return
     */
    public double detSlow() {
        if(this.m == 2) {
            return A[0][0]*A[1][1] - A[1][0]*A[0][1];
        } else {
            double sum = 0;
            double sign = 1;
            for(int j = 0; j < this.m; j++) {
                sum += sign*A[0][j]*detMatrix(j).detSlow();
                sign *= -1;
            }
            return sum;
        }
    }

    /**
     * Determinant
     * @return
     */
    public double det2() {
        if(this.m == 2) {
            return A[0][0]*A[1][1] - A[1][0]*A[0][1];
        } else {
            double sum = 0;
            double sign = 1;
            boolean[] all = new boolean[this.m]; 
            for(int i = 0; i < this.m; i++) {
                all[i] = true;
                sum += sign*A[0][i]*this.det2(1,all);
                sign *= -1;
                all[i] = false;
            }
            return sum;
        }
    }

    public double det2(int i, boolean[] skip) {
        if(this.m-i == 2) {
            int lo = 0;
            int hi = this.m;
            for(int j = lo; j < hi; j++) {
                if (!skip[j]) {
                    lo = j;
                    break;
                }
            }
            for(int j = lo+1; j < hi; j++) {
                if (!skip[j]) {
                    hi = j;
                    break;
                }
            }
            return A[i][lo]*A[i+1][hi] - A[i+1][lo]*A[i][hi];
        } else {
            double sum = 0;
            double sign = 1;
            for(int j = 0; j < this.m; j++) {
                if (!skip[j]) {
                    skip[j] = true;
                    sum += sign*A[i][j]*this.det2(i+1,skip);
                    sign *= -1;
                    skip[j] = false;
                }
            }
            return sum;
        }
    }

    public SquareMatrix LUMatrix() {
        SquareMatrix LU = new SquareMatrix(this.m);

        return LU;
    }

    public Object[] PLUDecompose() throws MatrixException {
        /*
        UpperTriangularMatrix U = new UpperTriangularMatrix(this.m);
        LowerTriangularMatrix L = new LowerTriangularMatrix(this.m);
        */
        int i,j,k,imax;

        double[][] A = new double[this.m][this.m];
        for(i = 0; i < this.m;i++) {
            System.arraycopy(this.A[i], 0, A[i], 0, this.m);
        }

        int[] P = new int[this.m+1];

        for (i = 0; i <= this.m; i++)
            P[i] = i;
            
        double maxA,absA;
        for (i = 0; i < this.m; i++) {
            maxA = 0.0;
            imax = i;
    
            for (k = i; k < this.m; k++)
                if ((absA = Math.abs(A[k][i])) > maxA) { 
                    maxA = absA;
                    imax = k;
                }
    
            if (maxA < 1e-120) {
                throw new MatrixException("Matrix is degenerate");
            } 
    
            if (imax != i) {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;
    
                //pivoting rows of A
                double[] row = A[i];
                A[i] = A[imax];
                A[imax] = row;
    
                //counting pivots starting from N (for determinant)
                P[this.m]++;
            }
    
            for (j = i + 1; j < this.m; j++) {
                A[j][i] /= A[i][i];
    
                for (k = i + 1; k < this.m; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
        }

        return new Object[]{P,new SquareMatrix(A)};
    }

    public double det() {
        if(this.m == 2) {
            return A[0][0]*A[1][1] - A[1][0]*A[0][1];
        } else {
            double sum = 0;
            double sign = 1;
            boolean[] all = new boolean[this.m]; 
            for(int i = 0; i < this.m; i++) {
                all[i] = true;
                sum += sign*A[0][i]*this.det2(1,all);
                sign *= -1;
                all[i] = false;
            }
            return sum;
        }
    }

    /**
     * Matrix minor (values exluding row i and col j)
     * @param j
     * @return
     */
    public SquareMatrix detMatrix(int j) {
        SquareMatrix out = new SquareMatrix(this.m-1);
        int col;
        for(int row = 1; row < this.m; row++) {
            int i = row-1;
            col = 0;
            while(col < j) {
                out.A[i][col] = this.A[row][col];
                col++;
            }
            col++;
            while(col < this.m) {
                out.A[i][col-1] = this.A[row][col];
                col++;
            }
        }
        return out;
    }

    /**
     * Matrix minor (determinant of submatrix exluding row i and col j)
     * @param j
     * @return
     */
    public double minor(int i, int j) {
        SquareMatrix out = new SquareMatrix(this.m-1);
        int row = 0;
        while(row < i) {
            for(int col = 0; col < m; col++) {
                if( col < j) {
                    out.A[row][col] = this.A[row][col];
                } else if(col > j) {
                    out.A[row][col-1] = this.A[row][col];
                }
            }
            row++;
        }
        row++;
        while(row < m) {
            for(int col = 0; col < m; col++) {
                if(col < j) {
                    out.A[i][col] = this.A[row][col];
                } else if(col > j) {
                    out.A[i][col-1] = this.A[row][col];
                }
            }
            i++;
            row++;
        }
        return out.det();
    }

    /**
     * Inverse of matrix
     * @return
     */
    public SquareMatrix inv() {
        double det = 0;
        SquareMatrix C = this.cofactor();
        for (int j = 0; j < this.m; j++) {
            det += this.A[0][j]*C.A[0][j];
        }
        det = 1/det;
        return new SquareMatrix(C.transpose().mult(det));
    }

    /**
     * Creates Matrix of minors
     * @return
     */
    public SquareMatrix cofactor() {
        SquareMatrix out = new SquareMatrix(this.m);
        double row_sign = 1;
        for(int i = 0; i < this.m; i++) {
            double sign = row_sign;
            for(int j = 0; j < this.m; j++) {
                out.A[i][j] = sign*this.minor(i,j);
                sign *= -1;
            }
            row_sign *= -1;
        }
        return out;
    }

    public Matrix LUSolve(Matrix b) throws MatrixException {
        double[] x = new double[b.m];

        int[] P = new int[this.m+1];

        int i,j,k,imax;

        for (i = 0; i <= this.m; i++)
            P[i] = i;
            
        double maxA,absA;
        for (i = 0; i < this.m; i++) {
            maxA = 0.0;
            imax = i;
    
            for (k = i; k < this.m; k++)
                if ((absA = Math.abs(A[k][i])) > maxA) { 
                    maxA = absA;
                    imax = k;
                }
    
            if (maxA < 1e-120) {
                throw new MatrixException("Matrix is degenerate");
            } 
    
            if (imax != i) {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;
    
                //pivoting rows of A
                double[] row = A[i];
                A[i] = A[imax];
                A[imax] = row;
    
                //counting pivots starting from N (for determinant)
                P[this.m]++;
            }
    
            for (j = i + 1; j < this.m; j++) {
                A[j][i] /= A[i][i];
    
                for (k = i + 1; k < this.m; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
        }

        for (i = 0; i < this.m; i++) {
            x[i] = b.A[P[i]][0];
    
            for (k = 0; k < i; k++)
                x[i] -= A[i][k] * x[k];
        }
    
        for (i = this.m - 1; i >= 0; i--) {
            for (k = i + 1; k < this.m; k++)
                x[i] -= A[i][k] * x[k];
    
            x[i] = x[i] / A[i][i];
        }

        return new Matrix(x,true);
    }

    public static Matrix GaussSeidel(SquareMatrix A, Matrix b, double tol, Matrix x) {
        double err = tol*2;
        int iter = 1;
        while(err > tol && iter < 10000)  {
            err = 0;
            for(int i = 0; i < b.m; i++) {
                double sum = b.A[i][0];
                for(int j = 0; j < b.m; j++) {
                    if(j != i) {
                        sum -= A.A[i][j]*x.A[j][0];
                    }
                }
                double tmp = x.A[i][0];
                x.A[i][0] = sum/A.A[i][i];
                tmp -= x.A[i][0];
                err += tmp*tmp;
            }
            iter++;
        }
        System.out.println("Solution in " +iter);
        return x;
    }

    public static Matrix SOR(SquareMatrix A, Matrix b, double tol, Matrix x, double relaxation) {
        double err = tol*2;
        int iter = 1;
        double w = 1-relaxation;
        while(err > tol && iter < 100000)  {
            err = 0;
            for(int i = 0; i < b.m; i++) {
                double sum = b.A[i][0];
                for(int j = 0; j < b.m; j++) {
                    if(j != i) {
                        sum -= A.A[i][j]*x.A[j][0];
                    }
                }
                double tmp = x.A[i][0];
                x.A[i][0] = w*x.A[i][0] + relaxation*sum/A.A[i][i];
                tmp -= x.A[i][0];
                err += tmp*tmp;
            }
            iter++;
        }
        System.out.println("Solution in " +iter);
        return x;
    }

    public static Matrix SOR(SquareMatrix A, Matrix b, double tol, Matrix x, double relaxationMin, double relaxationMax) {
        double err = tol*2;
        int iter = 1;
        double dw = (relaxationMax-relaxationMin)/100000;
        while(err > tol && iter < 100000)  {
            err = 0;
            double w = relaxationMin + dw*iter; 
            double w1 = 1-w;
            for(int i = 0; i < b.m; i++) {
                double sum = b.A[i][0];
                for(int j = 0; j < b.m; j++) {
                    if(j != i) {
                        sum -= A.A[i][j]*x.A[j][0];
                    }
                }
                double tmp = x.A[i][0];
                x.A[i][0] = w1*x.A[i][0] + w*sum/A.A[i][i];
                tmp -= x.A[i][0];
                err += tmp*tmp;
            }
            iter++;
        }
        System.out.println("Solution in " +iter);
        return x;
    }
}