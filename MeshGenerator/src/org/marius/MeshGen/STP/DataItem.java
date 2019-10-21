package org.marius.MeshGen.STP;

import java.util.ArrayList;

public class DataItem {
	public int id;
	public String type;
	public String args;
	
	public String name;
	
	ArrayList<Param> params = new ArrayList<Param>();
	
	public DataItem(){}

	public void parseArgs(){
		String arg = "";
		int i = 0;
		while(i < args.length()){
			int type = 0;
			char c = args.charAt(i);
			char flag = ',';
			if(c == '\'') {
				type = 1;
			} else if(c == '#') {
				type = 2;
			} else if(c == '.') {
				type = 3;
			} else if(c == '(') {
				type = 4;
				flag = ')';
				i++;
			} else if(c == '*') {
				type = 5;
				i++;
			} else {
				
			}
			boolean reading = true;
			
			while(reading) {
			}
				
			}
			i++;
		}
	}
	
	public String getClassName() {
		String[] temp = type.split("_");
		String classname = "";
		for(int i = 0; i < temp.length; i++) {
			classname += temp[i].charAt(0) + temp[i].substring(1).toLowerCase();
		}
		return classname;
	}
	
	public DataItem createSpecific() {
		try {
			Class<?> cls = Class.forName("org.marius.MeshGen.STP."+getClassName());
			try {
				DataItem entity = (DataItem) cls.newInstance();
				return entity;
			} catch (InstantiationException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		} catch (ClassNotFoundException e) {
			System.out.println("error matching entity. Generic Data Entity");
		}
		return new DataItem();
	}
}
