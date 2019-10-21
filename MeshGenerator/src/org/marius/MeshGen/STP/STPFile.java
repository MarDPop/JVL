package org.marius.MeshGen.STP;

import java.util.ArrayList;

import org.marius.MeshGen.Geometry.*;

public class STPFile {
	
	public static final String[] validGeometry = new String[]{""};

	private String name;
	
	public ArrayList<DataItem> dataEntities = new ArrayList<DataItem>();
	
	public STPFile(){};
	
	public STPFile(String filename) {
		this.name = filename;
	}
	
	public Geometry createGeometryFile() {
		Geometry g = new Geometry();
		
		return g;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public String getName() {
		return name;
	}
}
