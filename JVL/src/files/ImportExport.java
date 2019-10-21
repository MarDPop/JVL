package files;

import java.util.ArrayList;
import java.io.*;
import geometry.*;

public class ImportExport {

    public static void writeToXML(ArrayList<Surface> surfaces, String filename) {
        try {
            // Assume default encoding.
            FileWriter fileWriter =
                new FileWriter(filename+".xml");

            // Always wrap FileWriter in BufferedWriter.
            BufferedWriter bufferedWriter =
                new BufferedWriter(fileWriter);

            // Note that write() does not automatically
            // append a newline character.
            bufferedWriter.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            bufferedWriter.newLine();
            for(Surface s : surfaces) {
                bufferedWriter.write("<Surface>");
                bufferedWriter.newLine();
                int count = 0;
                for(Panel[] row : s.getPanels()) {
                    bufferedWriter.write("<Span id=\""+ count++ +"\">");
                    bufferedWriter.newLine();
                    for(Panel p : row) {
                        bufferedWriter.write("<Panel>");
                        bufferedWriter.newLine();
                        
                        bufferedWriter.write("</Panel>");
                        bufferedWriter.newLine();
                    }
                    bufferedWriter.write("</Span>");
                    bufferedWriter.newLine();
                }
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
    
}