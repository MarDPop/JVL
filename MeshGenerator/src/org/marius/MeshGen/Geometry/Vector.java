package org.marius.MeshGen.Geometry;

public class Vector {

	public final double[] x = new double[3];
	
	public Vector(){};
	
	public Vector(double[] y) {
		System.arraycopy(y, 0, this.x, 0, 3);
	}
	
	public Vector(Vector v) {
		System.arraycopy(v.x, 0, this.x, 0, 3);
	}
	
	public Vector(double x, double y, double z) {
		this.x[0] = x;
		this.x[1] = y;
		this.x[2] = z;
	}
	
	public void zero() {
		this.x[0] = 0;
		this.x[1] = 0;
		this.x[2] = 0;
	}
	
	public Vector add(Vector v) {
		Vector out = new Vector();
		out.x[0] = this.x[0] + v.x[0];
		out.x[1] = this.x[1] + v.x[1];
		out.x[2] = this.x[2] + v.x[2];
		return out;
	}
	
	public Vector subtract(Vector v) {
		Vector out = new Vector();
		out.x[0] = this.x[0] - v.x[0];
		out.x[1] = this.x[1] - v.x[1];
		out.x[2] = this.x[2] - v.x[2];
		return out;
	}
	
	public Vector add(double a) {
		Vector out = new Vector();
		out.x[0] = this.x[0] + a;
		out.x[1] = this.x[1] + a;
		out.x[2] = this.x[2] + a;
		return out;
	}
	
	public Vector multiply(double a) {
		Vector out = new Vector();
		out.x[0] = this.x[0]*a;
		out.x[1] = this.x[1]*a;
		out.x[2] = this.x[2]*a;
		return out;
	}
	
	public void addTo(Vector v) {
		this.x[0] += v.x[0];
		this.x[1] += v.x[1];
		this.x[2] += v.x[2];
	}
	
	public void multiplyBy(double a) {
		this.x[0] *= a;
		this.x[1] *= a;
		this.x[2] *= a;
	}
	
	public double dot(Vector v) {
		return this.x[0]*v.x[0]+this.x[1]*v.x[1]+this.x[2]*v.x[2];
	}
	
	public Vector cross(Vector v) {
		Vector out = new Vector();
		out.x[0] = this.x[1]*v.x[2] - this.x[2]*v.x[1];
		out.x[1] = this.x[2]*v.x[0] - this.x[0]*v.x[2];
		out.x[2] = this.x[0]*v.x[1] - this.x[1]*v.x[0];
		return out;
	}
	
	public double norm() {
		return Math.sqrt(this.dot(this));
	}
	
}
