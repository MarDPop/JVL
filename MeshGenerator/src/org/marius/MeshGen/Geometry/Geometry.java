package org.marius.MeshGen.Geometry;

import java.util.ArrayList;

public class Geometry {
	
	public String name;
	
	public ArrayList<Surface> surfaces = new ArrayList<Surface>();
	
	public Geometry(){};
	
	public void addSurface(Surface s) {
		surfaces.add(s);
	}
	
	public ArrayList<Surface> getSurfaces() {
		return surfaces;
	}
	
	public void toXML(String filename) {
		
	}
}
