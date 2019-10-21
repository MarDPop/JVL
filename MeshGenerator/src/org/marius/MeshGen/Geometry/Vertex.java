package org.marius.MeshGen.Geometry;

import java.util.HashSet;

public class Vertex {

	public Vector location;
	
	public HashSet<Face> faceList = new HashSet<Face>();
	
	public Vertex(){}
	
	public Vertex(double x, double y, double z) {
		location = new Vector(x,y,z);
	}
	
	public Vertex(double[] c) {
		location = new Vector(c);
	}
	
}
