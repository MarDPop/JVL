package org.marius.MeshGen;

import org.marius.MeshGen.Geometry.Geometry;
import org.marius.MeshGen.STP.STPFile;

public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	
		ImportTool importer = new ImportTool();
		System.out.println("working");
		STPFile stp = importer.importStep("/home/e385957/Desktop/basic.stp");
		Geometry g = stp.createGeometryFile();
		System.out.println("done");
	}

}
