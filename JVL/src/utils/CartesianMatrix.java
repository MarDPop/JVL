package utils;

/**
 * this matrix is the 3x3 
 */

 public class CartesianMatrix extends SquareMatrix {

    public CartesianMatrix() {
        super(3);
    }

    public CartesianMatrix(CartesianMatrix that) {
        super(3);
        for (int i = 0; i < 3; i++) {
            System.arraycopy(that.A[i], 0, this.A[i], 0, 3);
        }
    }

    public CartesianMatrix(Matrix B) {
        super(3);
        for(int i = 0; i < 3; i++) {
            System.arraycopy(B.A[i], 0, this.A[i], 0, 3);
        }
    }

    
    public CartesianMatrix(double yaw, double pitch, double roll) {
        super(3);
        rotationMatrixIntrinsic(yaw,pitch,roll);
    }
    
    public static CartesianMatrix getYXRotation(double yaw, double pitch) {
        CartesianMatrix out = new CartesianMatrix();
        double c1 = Math.cos(yaw);
        double s1 = Math.sin(yaw);
        double c2 = Math.cos(pitch);
        double s2 = Math.sin(pitch);
        out.A[0][0] = c1;
        out.A[0][1] = 0;
        out.A[0][2] = s1;
        out.A[1][0] = s1*s2;
        out.A[1][1] = c2;
        out.A[1][2] = -c1*s2;
        out.A[2][0] = -c2*s1;
        out.A[2][1] = s2;
        out.A[2][2] = c1*c2;
        return out;
    }

    public static CartesianMatrix getXYRotation(double yaw, double pitch) {
        CartesianMatrix out = new CartesianMatrix();
        double c1 = Math.cos(yaw);
        double s1 = Math.sin(yaw);
        double c2 = Math.cos(pitch);
        double s2 = Math.sin(pitch);
        out.A[0][0] = c1;
        out.A[0][1] = s1*s2;
        out.A[0][2] = s1*c2;
        out.A[1][0] = 0;
        out.A[1][1] = c2;
        out.A[1][2] = -s2;
        out.A[2][0] = -s1;
        out.A[2][1] = c1*s2;
        out.A[2][2] = c1*c2;
        return out;
    }

    public void getXRotationMatrix(double angle) {
        A[0][0] = 1;
        A[1][1] = Math.cos(angle);
        A[1][2] = -Math.sin(angle);
        A[2][1] = -A[2][1];
        A[2][2] = A[1][1];
    }

    public void getYRotationMatrix(double angle) {
        A[1][1] = 1;
        A[0][0] = Math.cos(angle);
        A[0][2] = Math.sin(angle);
        A[2][0] = -A[2][1];
        A[2][2] = A[1][1];
    }

    public void getZRotationMatrix(double angle) {
        A[2][2] = 1;
        A[0][0] = Math.cos(angle);
        A[0][1] = -Math.sin(angle);
        A[1][0] = -A[2][1];
        A[1][1] = A[1][1];
    }

    public void rotationMatrixIntrinsic(double yaw, double pitch, double roll) {
        double c1 = Math.cos(yaw);
        double s1 = Math.sin(yaw);
        double c2 = Math.cos(pitch);
        double s2 = Math.sin(pitch);
        double c3 = Math.cos(roll);
        double s3 = Math.sin(roll);
        A[0][0] = c1*c2;
        A[0][1] = c1*s2*s3-c3*s1;
        A[0][2] = s1*s3+c1*c3*s2;
        A[1][0] = c2*s1;
        A[1][1] = c1*c3+s1*s2*s3;
        A[1][2] = c3*s1*s2-c1*s3;
        A[2][0] = -s2;
        A[2][1] = c2*s3;
        A[2][2] = c2*c3;
    }

    @Override
    public CartesianMatrix inv() {
        CartesianMatrix out = new CartesianMatrix();
        out.A[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1];
        out.A[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0]);
        out.A[2][0] = A[1][0]*A[2][1]-A[1][1]*A[2][0];
        double det = 1/(A[0][0]*out.A[0][0]+A[0][1]*out.A[1][0]+A[0][2]*out.A[2][0]);
        out.A[0][0] *= det;
        out.A[1][0] *= det;
        out.A[2][0] *= det;

        out.A[0][1] = -det*(A[0][1]*A[2][2]-A[0][2]*A[2][1]);
        out.A[1][1] = det*(A[0][0]*A[2][2]-A[0][2]*A[2][0]);
        out.A[2][1] = -det*(A[0][0]*A[2][1]-A[0][1]*A[2][0]);

        out.A[0][2] = det*(A[0][1]*A[1][2]-A[0][2]*A[1][1]);
        out.A[1][2] = -det*(A[0][0]*A[1][2]-A[0][2]*A[1][0]);
        out.A[2][2] = det*(A[0][0]*A[1][1]-A[0][1]*A[1][0]);

        return out;
    }

    public Cartesian mult(Cartesian x) {
        Cartesian b = new Cartesian();
        b.x = A[0][0]*x.x + A[0][1]*x.y +A[0][2]*x.z ;
        b.y = A[1][0]*x.x + A[1][1]*x.y +A[1][2]*x.z ;
        b.z = A[2][0]*x.x + A[2][1]*x.y +A[2][2]*x.z ;
        return b;
    }
 }