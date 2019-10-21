package utils;

public final class MyMath {

    public static final double RAD2DEG = 180.0/Math.PI;
    public static final double DEG2RAD = Math.PI/180.0;
    public static final double GAS_R = 8.314459;
    public static final double HALFPI = Math.PI/2;
    public static final double TWOPI = Math.PI*2;
    public static final double TWOPI_INV = 0.159154943091895;

    // http://www.java-gaming.org/index.php?topic=36467.0
    private static final int Size_Ac = 100000;
    private static final int Size_Ar = Size_Ac + 1;
    private static final float Pi = (float) Math.PI;
    private static final float Pi_H = Pi / 2;

    private static final float Atan2[] = new float[Size_Ar];
    private static final float Atan2_PM[] = new float[Size_Ar];
    private static final float Atan2_MP[] = new float[Size_Ar];
    private static final float Atan2_MM[] = new float[Size_Ar];

    private static final float Atan2_R[] = new float[Size_Ar];
    private static final float Atan2_RPM[] = new float[Size_Ar];
    private static final float Atan2_RMP[] = new float[Size_Ar];
    private static final float Atan2_RMM[] = new float[Size_Ar];

    static {
        for (int i = 0; i <= Size_Ac; i++) {
            double d = (double) i / Size_Ac;
            double x = 1;
            double y = x * d;
            float v = (float) Math.atan2(y, x);
            Atan2[i] = v;
            Atan2_PM[i] = Pi - v;
            Atan2_MP[i] = -v;
            Atan2_MM[i] = -Pi + v;

            Atan2_R[i] = Pi_H - v;
            Atan2_RPM[i] = Pi_H + v;
            Atan2_RMP[i] = -Pi_H + v;
            Atan2_RMM[i] = -Pi_H - v;
        }
    }

    public static final float atan2(float y, float x) {
        if (y < 0) {
            if (x < 0) {
                //(y < x) because == (-y > -x)
                if (y < x) {
                    return Atan2_RMM[(int) (x / y * Size_Ac)];
                } else {
                    return Atan2_MM[(int) (y / x * Size_Ac)];
                }
            } else {
                y = -y;
                if (y > x) {
                    return Atan2_RMP[(int) (x / y * Size_Ac)];
                } else {
                    return Atan2_MP[(int) (y / x * Size_Ac)];
                }
            }
        } else {
            if (x < 0) {
                x = -x;
                if (y > x) {
                    return Atan2_RPM[(int) (x / y * Size_Ac)];
                } else {
                    return Atan2_PM[(int) (y / x * Size_Ac)];
                }
            } else {
                if (y > x) {
                    return Atan2_R[(int) (x / y * Size_Ac)];
                } else {
                    return Atan2[(int) (y / x * Size_Ac)];
                }
            }
        }
    }

    //https://en.wikipedia.org/wiki/Bhaskara_I%27s_sine_approximation_formula

    /**
     * faster cosine max err 9.1e-4
     * @param x
     * @return
     */
    public static double cos(double x) {
        if (Math.abs(x) > 6.283185307179586) {
            return qcos(x % 6.283185307179586);
        } else {
            if(Math.abs(x) < 1.570796326794897) {
                x*=x;
                return 1-x*(0.496072017833542-x*0.0367946); //0.999405-pi*pi*0.495589 + pi^4*0.036795
            } else {
                if(Math.abs(x) > 3.141592653589793) {
                    return -qcos(x-3.141592653589793);
                } else {
                    return -qcos(3.141592653589793-x);
                }
            }
        }
    }

    /**
     * faster cosine max err 5.2e-3 assumes x < +/- pi
     * @param x
     * @return
     */
    public static double qcos(double x) {
        x = Math.abs(x);
        if(x < 1.570796326794897) {
            return (1.570796326794897-x)/(1.570796326794897-x/(0.95013 - x));
        } else {
            return (1.570796326794897-x)*(4.091722653589793-x)/(3.285670260932521-0.570796326794897*x);
        }
    }

    public static double acos(double x) {
        boolean negate = x < 0;
        if (negate) {
            x = -x;
        }
        double ret = -0.0187293;
        ret = ret * x;
        ret = ret + 0.0742610;
        ret = ret * x;
        ret = ret - 0.2121144;
        ret = ret * x;
        ret = ret + 1.5707288;
        ret = ret * Math.sqrt(1.0-x);
        if (negate) {
            return ret - 3.14159265358979;
        } else {
            return 3.14159265358979 - ret;
        }
    }

    // http://www.coranac.com/2009/07/sines/

    public static double asin2(double x) {
        double x2 = x*x;
        return x*(1+x2*(0.166666666 + x2*(0.075 + x2*0.044642857142857)));
    }

    public static double asin(double x) {
        return x*(0.954929658551372 - 0.129006137732798*x*x);
    }

}