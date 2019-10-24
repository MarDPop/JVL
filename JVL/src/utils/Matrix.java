package utils;

public class Matrix {

    /**
     * values
     */
    protected final double[][] A;

    /**
     * number of rows
     */
    protected int m;

    /**
     * number of columns
     */
    protected int n;

    /**
     * May throw error if B rows are not same length
     * @param B
     */
    public Matrix(double[][] B) {
        this.A = B;
        this.m = B.length;
        this.n = B[0].length;
    }

    public Matrix(double[] B, boolean column) {
        if(column) {
            this.m = B.length;
            this.n = 1;
        } else {
            this.n = B.length;
            this.m = 1;
        }
        this.A = new double[m][n];
        for(int i = 0; i < B.length; i++) {
            if(column) {
                A[i][0] = B[i];
            } else {
                A[0][i] = B[i];
            }
        }
    }

    public Matrix(Matrix B) {
        this.A = new double[B.m][B.n];
        set(B);
    }

    public Matrix(int m, int n) {
        this.m = m;
        this.n = n;
        this.A = new double[m][n];
    }

    public void set(Matrix B) {
        this.m = B.A.length;
        this.n = B.A[0].length;
        for(int i = 0; i < m;i++) {
            System.arraycopy(B.A[i], 0, this.A[i], 0, n);
        }
    }

    /**
     * Gets value at index i,j
     * @param i
     * @param j
     * @return 
     */
    public double get(int i, int j) {
        return A[i][j];
    }

    /**
     * Sets value at index i,j
     * @param i
     * @param j
     * @param v
     */
    public void set(int i, int j, double v) {
        this.A[i][j] = v;
    }

    /**
     * Addition
     * @param that
     * @return
     */
    public Matrix add(Matrix that) {
        Matrix out = new Matrix(m,n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                out.A[i][j] = this.A[i][j] + that.A[i][j];
            }
        }
        return out;
    }

    /**
     * Subtraction
     * @param that
     * @return
     */
    public Matrix sub(Matrix that) {
        Matrix out = new Matrix(m,n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                out.A[i][j] = this.A[i][j] - that.A[i][j];
            }
        }
        return out;
    }

    /**
     * Matrix multiplication
     * @param that
     * @return
     */
    public Matrix mult(Matrix that) {
        Matrix out = new Matrix(this.n,that.m);
        for (int i = 0; i < this.m; i++) {
            for (int j = 0; j < that.n; j++) {
                for (int k = 0; k < that.m; k++) {
                    out.A[i][j]  += this.A[i][k]*that.A[k][j];
                }
            }
        }
        return out;
    }

    /** 
     * Scalar multiplication
     */
    public Matrix mult(double b) {
        Matrix out = new Matrix(m,n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                out.A[i][j] = this.A[i][j] *b;
            }
        }
        return out;
    }

    /**
     * Gets matrix values from source i1,j1 to i2,j2 exclusive
     * @param i1
     * @param j1
     * @param i2
     * @param j2
     * @return
     */
    public Matrix subMatrix(int i1, int j1, int i2, int j2) {
        Matrix out = new Matrix(i2-i1,j2-j1);
        for(int i = 0; i < out.m; i++) {
            int row = i1+i;
            for(int j = 0; j < out.n; j++) {
                out.A[i][j] = this.A[row][j1+j];
            }
        }
        return out;
    }

    /**
     * Matrix minor (values exluding row i and col j)
     * @param i
     * @param j
     * @return
     */
    public Matrix subMatrix(int i,int j) {
        Matrix out = new Matrix(this.m-1,this.n-1);
        int row = 0;
        while(row < i) {
            for(int col = 0; col < n; col++) {
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
            for(int col = 0; col < n; col++) {
                if( col < j) {
                    out.A[i][col] = this.A[row][col];
                } else if(col > j) {
                    out.A[i][col-1] = this.A[row][col];
                }
            }
            i++;
            row++;
        }
        return out;
    }

    public Matrix transpose() {
        Matrix out = new Matrix(this.n,this.m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                out.A[j][i] = this.A[i][j];
            }
        }
        return out;
    }
    
    public double[][] data() {
        return A;
    }

    /**
     * Gram Schmidt orthogonalization
     */
    public static Matrix gramSchmidt(Matrix V) {
        Matrix U = new Matrix(V);
        double sum_den, sum_num;
        for(int i = 0; i < V.n; i++) {

            for(int j = 0; j < i; j++) {
                sum_den = 0;
                sum_num = 0;
                for(int k = 0; k < U.m; j++) {
                    sum_den += U.A[k][j]*U.A[k][j];
                    sum_num += U.A[k][i]*U.A[k][j];
                }
                sum_num /= sum_den;
                for(int k = 0; k < U.m; j++) {
                    U.A[k][i] -= sum_num*U.A[k][j];
                }
            }

            sum_den = 0;
            for(int j = 0; j < U.m; j++) {
                sum_den += U.A[j][i]*U.A[j][i];
            }
            sum_den = 1/Math.sqrt(sum_den);
            for(int j = 0; j < U.m; j++) {
                U.A[j][i] *= sum_den;
            }
        }
        return U;
    }

    /**
     * QR decomposition
     * @return
     */
    public static Matrix[] GivensQR(Matrix A) {
        double a,b;
        SquareMatrix Q = new SquareMatrix(A.m,1);
        Matrix R = new Matrix(A);
        for (int j = 0; j < A.n; j++) {
            for (int i = A.m; i > j; i--) {
                double[] coef = MyMath.Givins(R.A[i-1][j],R.A[i][j]);
                R.A[i-1][j] = coef[2]; 

                for(int k = 0; k < A.m; k++){
                    a = Q.A[k][i-1];
                    b = Q.A[k][i];
                    Q.A[k][i-1] = a*coef[0] + coef[1]*b; 
                    Q.A[k][i] = coef[0]*b - coef[1]*a;  
                }
                
                for(int k = j+1; k < A.n; k++) {
                    a = R.A[i-1][k];
                    b = R.A[i][k];
                    R.A[i-1][k] = coef[0]*a + coef[1]*b; 
                    R.A[i][k] = coef[0]*b - coef[1]*a; 
                    
                }
            }
        }
        return new Matrix[]{Q,R};
    }

    /**
     * QR decomposition
     * @return
     */
    public static Matrix[] ParallelGivensQR(Matrix A) {
        double a,b;
        SquareMatrix Q = new SquareMatrix(A.m,1);
        Matrix R = new Matrix(A);
        for (int j = 0; j < A.n; j++) {
            for (int i = A.m; i > j; i--) {
                double[] coef = MyMath.Givins(R.A[i-1][j],R.A[i][j]);
                R.A[i-1][j] = coef[2]; 

                for(int k = 0; k < A.m; k++){
                    a = Q.A[k][i-1];
                    b = Q.A[k][i];
                    Q.A[k][i-1] = a*coef[0] + coef[1]*b; 
                    Q.A[k][i] = coef[0]*b - coef[1]*a;  
                }
                
                for(int k = j+1; k < A.n; k++) {
                    a = R.A[i-1][k];
                    b = R.A[i][k];
                    R.A[i-1][k] = coef[0]*a + coef[1]*b; 
                    R.A[i][k] = coef[0]*b - coef[1]*a; 
                    
                }
            }
        }
        return new Matrix[]{Q,R};
    }

    /**
     * QR decomposition
     * @return
     */
    public static Matrix[] HouseholderQR(Matrix A) {
        SquareMatrix Q = new SquareMatrix(A.m,1);
        Matrix R = new Matrix(A);
        double sum, u1, tau;
        double[] w = new double[A.m];
        double[][] T;
        for(int j = 0; j < A.n; j++) {
            sum = 0;
            for(int k = 0; k < A.m; k++){
                sum += R.A[k][j]*R.A[k][j];
            }
            tau = Math.signum(R.A[j][j])/Math.sqrt(sum);
            u1 = 1/(R.A[j][j] + 1/tau);
            w[j] = 1;
            for(int k = j+1; k < A.m;k++) {
                w[k] = R.A[k][j]*u1;
            }
            tau /= u1;

            T = new double[A.m-j][A.n];

            for(int k = j; k < A.m; k++) {
                for( int i = 0; i < A.n; i++) {
                    double r_dot = 0;
                    for (int t = j; t < A.m; t++){
                        r_dot += R.A[t][i]*w[t];
                    }
                    T[k-j][i] = tau*w[k]*r_dot;
                } 
            }

            for(int k = j; k < A.m; k++) {
                for( int i = 0; i < A.n; i++) {
                    R.A[k][i] -= T[k-j][i];
                } 
            }

            T = new double[A.m][A.m-j];

            for(int k = j; k < A.m; k++) {
                for( int i = 0; i < A.n; i++) {
                    double r_dot = 0;
                    for (int t = j; t < A.m; t++){
                        r_dot += Q.A[i][t]*w[t];
                    }
                    T[i][k-j] = tau*w[k]*r_dot;
                } 
            }

            for(int k = j; k < A.m; k++) {
                for( int i = 0; i < A.m; i++) {
                    Q.A[i][k] -= T[i][k-j];
                } 
            }
        }

        return new Matrix[]{Q,R};
    }

    public static double dot(Matrix A, Matrix B) {
        if(A.n == 1) {
            if(B.m == A.m) {
                double sum = 0;
                for(int i = 0; i < A.m; i++) {
                    sum += A.A[i][0]*B.A[i][0];
                }
                return sum;
            } else {
                return Double.NaN;
            }
        } else {
            if(A.n == B.n) {
                double sum = 0;
                for(int i = 0; i < A.n; i++) {
                    sum += A.A[0][i]*B.A[0][1];
                }
                return sum;
            } else {
                return Double.NaN;
            }
        }
    }

}