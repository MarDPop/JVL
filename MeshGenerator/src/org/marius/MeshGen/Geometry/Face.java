package org.marius.MeshGen.Geometry;

public class Face {
	
	public Vector center = new Vector();
	
	private final Vertex[] points = new Vertex[3];
	
	private Vector normal;
	
	public Face(){}
	
	public Face(Vertex[] points){
		create(points);
	}
	
	public void create(Vertex[] points) {
		System.arraycopy(points, 0, this.points, 0, 3);
		center.zero();
		for(int i = 0; i < 3; i++) {
			center.addTo(points[i].location);
			points[i].faceList.add(this);
		}
		center.multiplyBy(1/3);
		Vector v1 = points[1].location.subtract(points[0].location);
		Vector v2 = points[2].location.subtract(points[0].location);
		normal = v1.cross(v2);
		normal.multiply(0.5); // normal is area size
	}
	
	public void flipNormal() {
		this.normal.multiply(-1);
	}

}
