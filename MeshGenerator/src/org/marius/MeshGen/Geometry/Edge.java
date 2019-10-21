package org.marius.MeshGen.Geometry;

public class Edge {
	
	private Vector vector;
	
	private double length;

	public final Vertex[] points = new Vertex[2];
	
	public Edge(){}
	
	public void setPoints(Vertex v1, Vertex v2) {
		points[0] = v1;
		points[1] = v2;
		vector = v2.location.subtract(v1.location);
		length = vector.norm();
	}
	
	public double getLength() {
		return length;
	}
	
}
