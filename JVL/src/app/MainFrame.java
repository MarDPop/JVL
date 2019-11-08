package app;

import java.awt.*;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.awt.event.ActionEvent;

import javax.swing.*;
import app.view.ViewingPane;

public class MainFrame extends JFrame {

    /**
     *
     */
    private static final long serialVersionUID = 1729233823098776690L;

    public ViewingPane vp = new ViewingPane();

    public static JMenuBar bar = new JMenuBar();

    public static JMenu fileMenu = new JMenu("File");

    public static JMenu geometryMenu = new JMenu("Geometry");

    public static JMenuItem addGeometry = new JMenuItem("Add");

    public static JMenu simMenu = new JMenu("Simulation");

    public static JMenu help = new JMenu("Help");

    public static JMenuItem about = new JMenuItem("About");

    public static JMenuItem readMe = new JMenuItem("Read Me");

    public MainFrame(String version) {
        this.setTitle("JVL");

        fileMenu.add(new JMenuItem("Open"));
        fileMenu.add(new JMenuItem("Save"));

        
        geometryMenu.add(new JMenuItem("Import"));
        addGeometry.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                // addGeometryFrame();
            }
        });

        geometryMenu.add(new JMenuItem("Add"));
        geometryMenu.add(new JMenuItem("Edit"));

        simMenu.add(new JMenuItem("Setup"));
        simMenu.add(new JMenuItem("Run"));
        simMenu.add(new JMenuItem("Options"));

        about.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JOptionPane.showMessageDialog(null, "<html><p>JVL Version " + version  + "</p><p>Author: Marius Popescu </p></html>", "About", JOptionPane.INFORMATION_MESSAGE);
            }
        });

        readMe.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                File file = new File("README.html");
                try {
                    Desktop.getDesktop().open(file);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
        });

        help.add(about);
        help.add(readMe);

        bar.add(fileMenu);
        bar.add(geometryMenu);
        bar.add(simMenu);
        bar.add(help);

        getContentPane().add(BorderLayout.NORTH, bar);

        getContentPane().add(vp);

        vp.getResult = true;

    } 


    public static void addGeometryFrame() {
        EventQueue.invokeLater(new Runnable()
        {

			@Override
			public void run() {
                JFrame popup = new JFrame();

                popup.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                popup.setSize(500,500);
                
                JTextField OriginX = new JTextField(); 
                JTextField OriginY = new JTextField(); 
                JTextField OriginZ = new JTextField(); 
                popup.getContentPane().add(OriginX);
                popup.getContentPane().add(OriginY);
                popup.getContentPane().add(OriginZ);

                popup.setVisible(true);
			}
            
        });
    }
}