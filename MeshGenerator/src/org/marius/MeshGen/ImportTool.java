package org.marius.MeshGen;

import org.marius.MeshGen.STP.*;

import java.io.*;
import java.util.Scanner;
import java.util.regex.*;

public class ImportTool {
	
	public ImportTool(){}
	
	/*
	 * Uses ISO 10303-21 STEP file for protocol AP203
	 * 
	 */
	public STPFile importStep(String filePath){
		STPFile s = new STPFile();
		try 
	    {
			Scanner scanner = new Scanner(new FileReader(filePath));
			scanner.useDelimiter(";");
	        String line = null;
	        while (scanner.hasNext()) 
	        {
	        	line = scanner.next();
	        	
	            if(line.indexOf("FILE_NAME") >= 0) {
	            	int idx1 = line.indexOf('\'');
	            	int idx2 = line.indexOf('\'',idx1);
	            	s.setName(line.substring(idx1,idx2));
	            }
	            if(line.indexOf("DATA") >= 0) {
	            	break;
	            }
	        }
	        
	        Pattern format = Pattern.compile(".*#(\\d+)=(\\w+)\\((.*)\\)");
	        while (scanner.hasNext()) 
	        {
	        	line = scanner.next();
	        	Matcher m = format.matcher(line);
	        	if(m.find()) {
	        		int id = Integer.parseInt(m.group(1));
	        		String entityType = m.group(2);
	        		String[] temp = entityType.split("_");
	        		String classname = "";
	        		for(int i = 0; i < temp.length; i++) {
	        			classname += temp[i].charAt(0) + temp[i].substring(1).toLowerCase();
	        		}
	        		
	        		DataItem entity = new DataItem();
	        		
	        		entity.id = id;
					entity.type = entityType;
					entity.args = m.group(3);
	        		s.dataEntities.set(id,entity);
	        		
	        		System.out.println(classname);
	        	} else {
	        		System.out.println("didnt find match");
	        	}
	        	
	        }
	        
	        scanner.close();
	        
	    }
	    catch (IOException e) 
	    {
	        e.printStackTrace();
	    }
		
		return s;
	}

}
