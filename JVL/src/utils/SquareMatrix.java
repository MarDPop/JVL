package utils;

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

    public SquareMatrix(int m, double a) {
        super(m,m);
        for(int i = 0; i < m ; i++) {
            A[i][i] = -a;
        }
    }

    /**
     * Decomposes this matrix into a Pivoting array and Square Matrix consisting of Lower and Upper triangular decomposition
     */
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

    /**
     * Returns determinant
     * @return
     */
    public double det() {
        if(this.m == 2) {
            return A[0][0]*A[1][1] - A[1][0]*A[0][1];
        } else if (this.m == 3){
            return A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
        } else if (this.m < 7) {
            // technically 5 is largest where n! is about same as 2/3n^3 + n^2 but pivoting can be slow and this is more stable
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
        } else {
            try {
                return LUPDeterminant();
            } catch (MatrixException e) {
                return -1;
            }
        }
    }

    /**
     * Determinant helper function
     * @param i
     * @param skip
     * @return
     */
    protected double det2(int i, boolean[] skip) {
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
        if (this.m < 5) {
            double det = 0;
            SquareMatrix C = this.cofactor();
            for (int j = 0; j < this.m; j++) {
                det += this.A[0][j]*C.A[0][j];
            }
            det = 1/det;
            return new SquareMatrix(C.transpose().mult(det));
        } else {
            try {
                return LUInvert();
            } catch (MatrixException e) {
                return null;
            }
        }
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

    /**
     * Solves Ax = B using LU decomposition overwriting Matrix A
     * @param b
     * @return
     * @throws MatrixException
     */
    public Matrix LUOverwriteSolve(Matrix b) throws MatrixException {
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

    /**
     * Solves Ax = b but A is is LU matrix P is pivoting array
     * @param A
     * @param P
     * @param b
     * @param x
     */
    public static void LUSolve(SquareMatrix LU, int[] P, Matrix b, Matrix x) {
        int i,k;
        for (i = 0; i < LU.m; i++) {
            x.A[i][0] = b.A[P[i]][0];
    
            for (k = 0; k < i; k++)
            x.A[i][0] -= LU.A[i][k] * x.A[k][0];
        }
    
        for (i = LU.m - 1; i >= 0; i--) {
            for (k = i + 1; k < LU.m; k++)
            x.A[i][0] -= LU.A[i][k] * x.A[k][0];
    
            x.A[i][0] /= LU.A[i][i];
        }
    }

    /**
     * Invert Matrix using LU
     * @return
     * @throws MatrixException
     */
    public SquareMatrix LUInvert() throws MatrixException {
        Object[] arr = this.PLUDecompose();
        int[] P = (int[]) arr[0];
        double[][] A = ((SquareMatrix) arr[1]).A;
        double[][] IA = new double[this.m][this.m];
        for (int j = 0; j < this.m; j++) {
            for (int i = 0; i < this.m; i++) {
                if (P[i] == j) 
                    IA[i][j] = 1.0;
    
                for (int k = 0; k < i; k++)
                    IA[i][j] -= A[i][k] * IA[k][j];
            }
    
            for (int i = this.m - 1; i >= 0; i--) {
                for (int k = i + 1; k < this.m; k++)
                    IA[i][j] -= A[i][k] * IA[k][j];
    
                IA[i][j] = IA[i][j] / A[i][i];
            }
        }
        return new SquareMatrix(IA);
    }

    /**
     * Determine determinant using LU decomposition
     * @return
     * @throws MatrixException
     */
    public double LUPDeterminant() throws MatrixException {
        Object[] arr = this.PLUDecompose();

        int[] P = (int[]) arr[0];
        double[][] A = ((SquareMatrix) arr[1]).A;

        double det = A[0][0];
        
        for (int i = 1; i < this.m; i++)
            det *= A[i][i];
    
        if ((P[this.m] - this.m) % 2 == 0)
            return det; 
        else
            return -det;
    }

    public void invert() {
        
    }

    /* STATIC METHODS */

    /**
     * Iterative Gauss - Seidel for solution
     * @param A
     * @param b
     * @param tol
     * @param x
     * @return
     */
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

    /**
     * Successive Over Relaxation Solve
     * @param A
     * @param b
     * @param tol
     * @param relaxation
     * @return
     */
    public static Matrix SOR(SquareMatrix A, Matrix b, double tol, double relaxation) {
        double[] x = new double[b.m];
        double err = tol*2;
        int iter = 1;
        double w = 1-relaxation;
        while(err > tol && iter < 100000)  {
            err = 0;
            for(int i = 0; i < b.m; i++) {
                double sum = b.A[i][0];
                for(int j = 0; j < b.m; j++) {
                    if(j != i) {
                        sum -= A.A[i][j]*x[j];
                    }
                }
                double tmp = x[i];
                x[i] = w*x[i] + relaxation*sum/A.A[i][i];
                tmp -= x[i];
                err += tmp*tmp;
            }
            iter++;
        }
        System.out.println("Solution in " +iter);
        return new Matrix(x,true);
    }

    /**
     * Successive Over Relaxation method with a linear relaxation parameter
     * @param A
     * @param b
     * @param tol
     * @param x
     * @param relaxationMin
     * @param relaxationMax
     */
    public static void SOR(SquareMatrix A, Matrix b, double tol, Matrix x, double relaxationMin, double relaxationMax) {
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
    }

    /**
     * Conjugate Gradient Solve
     * @param A
     * @param b
     * @param tol
     * @param x
     */
    public static void CG(SquareMatrix A, Matrix b, double tol, Matrix x) {
        Matrix r = b.sub(A.mult(x));

        int i;
        double err = 0;
        for(i = 0; i < A.m;i++) {
            err += r.A[i][0]*r.A[i][0];
        }

        Matrix r_next = new Matrix(r);
        Matrix p = new Matrix(r);

        int iter = 0;
        
        double alpha,beta;
        while(err > tol && iter < 10000)  {
            Matrix tmp = A.mult(p);

            alpha = dot(r,r)/dot(p,tmp);
            err = 0;
            for(i = 0; i < A.m;i++){
                x.A[i][0] += alpha*p.A[i][0];
                r_next.A[i][0] = r.A[i][0] - alpha*tmp.A[i][0];
                err += r_next.A[i][0]*r_next.A[i][0];
            }

            if(err < tol)
                break;

            beta = dot(r_next,r_next)/dot(r,r);
            for(i = 0; i < A.m;i++){
                p.A[i][0] = r_next.A[i][0]+beta*p.A[i][0];
                r.A[i][0] = r_next.A[i][0];
            }
            iter++;
        }
        System.out.println("Solution in " +iter);
    }

    /**
     * Invert using Strassen Algorithm
     * @param A
     * @return
     */
    public static SquareMatrix StrassenInvert(SquareMatrix A) {
        return A;
    }

    /**
     * Singular Value Decomposition
     * @param A
     * @return
     */
    public static Object[] SVD(SquareMatrix A) {
        return null;
    }    

    public static Matrix powereig(SquareMatrix A, Matrix guess) {
        for(int iter = 0; iter < 1000; iter++) {
            Matrix b = A.mult(guess);
            Matrix diff = b.sub(guess);

            double err = 0;
            double norm = 0 ;
            for (int i = 0; i < A.m; i++) {
                err += diff.A[i][0]*diff.A[i][0];
                norm += b.A[i][0]*b.A[i][0];
            }
            norm = 1/Math.sqrt(norm);
            guess.set(b.mult(norm));
            if(err*norm < 1e-12) {
                return b;
            }
        }
        return null;
    }

    public static double rayleighQuotient(SquareMatrix A, Matrix b) {
        double a = Matrix.dot(b,A.mult(b));
        return a /= Matrix.dot(b,b);
    }

    public static Matrix inveig(SquareMatrix A, Matrix guess, double shift) throws MatrixException {
        SquareMatrix I = new SquareMatrix(A);
        for(int i = 0; i < A.m;i++) {
            I.A[i][i] -= shift;
        }
        
        Object[] arr = I.PLUDecompose();
        SquareMatrix LU = (SquareMatrix)arr[1];
        int[] P = (int[]) arr[0];

        Matrix w = new Matrix(guess);
        for(int iter = 0; iter < 1000; iter++) {
            double r = Matrix.normalize(guess);
            SquareMatrix.LUSolve(LU, P, guess, w);
            Matrix diff = w.sub(guess);
            if(Matrix.norm(diff) < 0.0001*r){
                return w;
            }
        }
        return null;
    }

    public static double eigenvalueFromEigenvector(SquareMatrix A, Matrix eigenvector) {
        normalize(eigenvector);
        return norm(A.mult(eigenvector));
    }

    /**
     * Givens Matrix for i and j where j > i
     * @return
     */
    public static SquareMatrix GivensMatrix(int i, int j, double theta, int size) {
        SquareMatrix I = new SquareMatrix(size,1);
        I.A[i][i] = Math.cos(theta);
        I.A[j][j] = I.A[i][i];
        I.A[j][i] = Math.sin(theta);
        I.A[i][j] = -I.A[j][i];
        return null;
    }

    /**
     * Cholesky Decomposition (use only for symmetric)
     */
    public static LowerTriangularMatrix chol(SquareMatrix A) {
        LowerTriangularMatrix L = new LowerTriangularMatrix(A.m);
        double sum;
        for(int k = 0; k < A.m; k++) {
            for(int i = 0; i < k; i++) {
                sum = 0;
                for(int j = 0; j < i; j++) {
                    sum += L.A[i][j]*L.A[k][j];
                }
                L.A[k][i] = A.A[k][i] - sum/L.A[i][i];
            }
            sum = 0;
            for(int j = 0; j < k; j++) {
                sum += L.A[k][j]*L.A[k][j];
            }
            L.A[k][k] = Math.sqrt(A.A[k][k]-sum);
        }
        return L;
    }
}