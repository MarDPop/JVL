package utils;

/**
 * A cartesian based 
 */
public class Cartesian {

    /**
     * 3 dimensional values
     */
    public double x;

    public double y;

    public double z;

    /**
     * euler norm
     */
    public double r;

    /* Constructors */
    public Cartesian() { }

    public Cartesian(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    /**
     * Quick setter
     * @param x
     */
    public Cartesian(double[] x) {
        this.x = x[0];
        this.y = x[1];
        this.z = x[2];
    }

    public Cartesian(Cartesian a) {
        this.x = a.x;
        this.y = a.y;
        this.z = a.z;
        this.r = a.r;
    }

    /**
     * Setter for values
     * @param x
     * @param y
     * @param z
     */
    public void set(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
        getMagnitude();
    }

    /**
     * 
     * @param copy
     */
    public void set(Cartesian copy) {
        this.x = copy.x;
        this.y = copy.y;
        this.z = copy.z;
        this.r = copy.r;
    }

    /**
     * dot product
     * @param a input vector
     * @return dot product
     */
    public double dot(Cartesian a) {
        return a.x*this.x+a.y*this.y+a.z*this.z;
    }

    /**
     * element multiply
     * @param a
     * @return
     */
    public Cartesian scale(Cartesian a) {
        Cartesian out = new Cartesian();
        out.x = this.x*a.x;
        out.y = this.y*a.y;
        out.z = this.z*a.z;
        return out;
    }

    public Cartesian scaleBy(Cartesian a) {
        this.x*=a.x;
        this.y*=a.y;
        this.z*=a.z;
        return this;
    }

    public Cartesian add(Cartesian a) {
        return new Cartesian(this.x+a.x,this.y+a.y,this.z+a.z);
    }

    public Cartesian addTo(Cartesian a) {
        this.x+=a.x;
        this.y+=a.y;
        this.z+=a.z;
        return this;
    }

    public Cartesian sub(Cartesian a) {
        return new Cartesian(this.x-a.x,this.y-a.y,this.z-a.z);
    }

    public Cartesian subFrom(Cartesian a) {
        this.x-=a.x;
        this.y-=a.y;
        this.z-=a.z;
        return this;
    }

    public Cartesian mult(double b) {
        return new Cartesian(this.x*b, this.y*b,this.z*b);
    }

    public Cartesian multBy(double b) {
        this.x*=b;
        this.y*=b;
        this.z*=b;
        return this;
    }

    public Cartesian mult(CartesianMatrix A) {
        Cartesian out = new Cartesian();
        out.x = A.get(0,0)*this.x + A.get(1,0)*this.y + A.get(2,0)*this.z;
        out.y = A.get(0,1)*this.x + A.get(1,1)*this.y + A.get(2,1)*this.z;
        out.z = A.get(0,2)*this.x + A.get(1,2)*this.y + A.get(2,2)*this.z;
        return out;
    }

    public Cartesian normalize() {
        getMagnitude();
        if(r == 0) {
            return this;
        }
        this.x /= this.r;
        this.y /= this.r;
        this.z /= this.r;
        this.r = 1;
        return this;
    }

    /**
     * returns euler norm
     */
    public double getMagnitude() {
        return this.r = Math.sqrt(x*x + y*y + z*z);
    }

    /**
     * cross product u x v
     * @param a vector v
     * @return 
     */
    public Cartesian cross(Cartesian a) {
        Cartesian out =  new Cartesian();
        out.x = this.y*a.z - this.z*a.y;
        out.y = this.z*a.x - this.x*a.z;
        out.z = this.x*a.y - this.y*a.x;
        return out;
    }

}