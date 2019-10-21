package utils;

public class SquareMatrix extends Matrix {

    public SquareMatrix(int m) {
        super(m,m);
    }

    public SquareMatrix(SquareMatrix that) {
        super(that);
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
    public double det() {
        if(m == 2) {
            return A[0][0]*A[1][1] - A[1][0]*A[0][1];
        } else {
            double sum = 0;
            double sign = 1;
            for(int j = 0; j < n; j++) {
                sum += sign*this.detMatrix(j).det();
                sign *= -1;
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
        for(int row = 1; row < m; row++) {
            int i = row-1;
            col = 0;
            while(col < j) {
                out.A[i][col] = this.A[row][col];
                col++;
            }
            col++;
            while(col < n) {
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
        for (int j = 0; j < n; j++) {
            det += C.A[0][j];
        }
        det = 1/det;
        return (SquareMatrix)C.transpose().mult(det);
    }

    /**
     * Creates Matrix of minors
     * @return
     */
    public SquareMatrix cofactor() {
        SquareMatrix out = new SquareMatrix(this.m);
        double sign = 1;
        for(int i = 0; i < this.m; i++) {
            for(int j = 0; j < this.m; j++) {
                out.A[i][j] = sign*this.minor(i,j);
                sign *= -1;
            }
            //sign *= -1;
        }
        return out;
    }

    public Matrix LUSolve(Matrix b) {
        Matrix x = new Matrix(b.m,1);
        return x;
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
}