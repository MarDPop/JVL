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
        this.m = B.A.length;
        this.n = B.A[0].length;
        this.A = new double[m][n];
        for(int i = 0; i < m;i++) {
            System.arraycopy(B.A[i], 0, this.A[i], 0, n);
        }
    }

    public Matrix(int m, int n) {
        this.m = m;
        this.n = n;
        this.A = new double[m][n];
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
}