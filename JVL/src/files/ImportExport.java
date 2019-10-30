package files;

import java.util.ArrayList;
import java.util.Date;
import java.io.*;
import geometry.*;

public class ImportExport {

    public static void writeToXML(ArrayList<Surface> surfaces, String filename) {
        try {
            // Assume default encoding.
            FileWriter fileWriter =
                new FileWriter(filename+".surf");

            // Always wrap FileWriter in BufferedWriter.
            BufferedWriter bufferedWriter =
                new BufferedWriter(fileWriter);

            // Note that write() does not automatically
            // append a newline character.
            bufferedWriter.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            bufferedWriter.newLine();
            for(Surface s : surfaces) {
                bufferedWriter.write("<Surface >");
                bufferedWriter.newLine();
                int count = 0;
                for(Panel[] row : s.getPanels()) {
                    bufferedWriter.write("<Span id=\""+ count++ +"\">");
                    bufferedWriter.newLine();
                    for(Panel p : row) {
                        bufferedWriter.write("<Panel>");
                        bufferedWriter.newLine();

                        bufferedWriter.write("<Vertex inboard='true' upstream='true'>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("<x>"+p.vertices[0].x+"</x>"+"<y>"+p.vertices[0].y+"</y>"+"<z>"+p.vertices[0].z+"</z>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("</Vertex>");
                        bufferedWriter.newLine();

                        bufferedWriter.write("<Vertex inboard='true' upstream='false'>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("<x>"+p.vertices[1].x+"</x>"+"<y>"+p.vertices[1].y+"</y>"+"<z>"+p.vertices[1].z+"</z>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("</Vertex>");
                        bufferedWriter.newLine();

                        bufferedWriter.write("<Vertex inboard='false' upstream='false'>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("<x>"+p.vertices[2].x+"</x>"+"<y>"+p.vertices[2].y+"</y>"+"<z>"+p.vertices[2].z+"</z>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("</Vertex>");
                        bufferedWriter.newLine();

                        bufferedWriter.write("<Vertex inboard='false' upstream='true'>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("<x>"+p.vertices[3].x+"</x>"+"<y>"+p.vertices[3].y+"</y>"+"<z>"+p.vertices[3].z+"</z>");
                        bufferedWriter.newLine();
                        bufferedWriter.write("</Vertex>");
                        bufferedWriter.newLine();

                        bufferedWriter.write("</Panel>");
                        bufferedWriter.newLine();
                    }
                    bufferedWriter.write("</Span>");
                    bufferedWriter.newLine();
                }

                bufferedWriter.write("<Controls>");
                bufferedWriter.newLine();
                for(int[] control : s.getControls()) {
                    bufferedWriter.write("<Definition type='"+control[0]+"'>");
                    bufferedWriter.newLine();
                    bufferedWriter.write("<StartRow>"+control[1]+"</StartRow>");
                    bufferedWriter.newLine();
                    bufferedWriter.write("<EndRow>"+control[2]+"</EndRow>");
                    bufferedWriter.newLine();
                    bufferedWriter.write("<StartCol isRotation='true'>"+control[3]+"</StartCol>");
                    bufferedWriter.newLine();
                    bufferedWriter.write("<EndCol>"+control[4]+"</EndCol>");
                    bufferedWriter.newLine();
                    bufferedWriter.write("<Deflection>"+control[5]+"</Deflection>");
                    bufferedWriter.newLine();
                    bufferedWriter.write("</Definition>");
                    bufferedWriter.newLine();
                }
                bufferedWriter.write("</Controls>");
                bufferedWriter.newLine();

                bufferedWriter.write("</Surface>");
                bufferedWriter.newLine();
            }
            
            
            // Always close files.
            bufferedWriter.close();
        }
        catch(IOException ex) {
            ex.printStackTrace();
        }
    }

    public static void writeCoefficients(double[] coefficients, String[] names, double[] references, String filename) {
        try {
            // Assume default encoding.
            FileWriter fileWriter =
                new FileWriter(filename);

            // Always wrap FileWriter in BufferedWriter.
            BufferedWriter bufferedWriter =
                new BufferedWriter(fileWriter);

            bufferedWriter.write("Coefficients from Run performed at "+ new Date().toString() );
            bufferedWriter.newLine();

            for(int i = 0; i < coefficients.length;i++) {
                if ( i < names.length) {
                    bufferedWriter.write(names[i] + " = ");
                } else {
                    bufferedWriter.write("? = ");
                }
                bufferedWriter.write("" + coefficients[i]);
                bufferedWriter.newLine();
            }

            bufferedWriter.close();
        } 
        catch(IOException ex) {
            ex.printStackTrace();
        }
    }
    
}