package utils;

public class TridiagonalMatrix extends Matrix {

    double a; // lower diagonal

    double b; // primary diagonal

    double c; // upper diagonal

    public TridiagonalMatrix(double a, double b, double c, int n) {
        super(3,n);
        this.n = n;
        this.a = a;
        this.b = b;
        this.c = c;
    }
}