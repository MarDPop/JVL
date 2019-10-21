package org.marius.MeshGen.Geometry;

import java.util.ArrayList;

public class Cell {
	
	public ArrayList<Face> faces = new ArrayList<Face>();
	
	public Vector center = new Vector();
	
	public ArrayList<Cell> neighbors = new ArrayList<Cell>();
	
	public ArrayList<Double> dx = new ArrayList<Double>();
	
	public Properties prop;
	
	public Cell(){}
	
}
