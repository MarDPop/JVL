package app;

import javax.swing.*;
import java.awt.*;
import geometry.*;
import utils.Cartesian;
import utils.MyMath;
import utils.SquareMatrix;
import utils.Matrix;

public class App {
    static JFrame gui;

    static VL2 sim = new VL2();

    static ViewingPane vp = new ViewingPane();

    public static void main(String[] args) throws Exception {
        simTest();
    }

    private static void simTest() {
        gui = new JFrame("JVL");
        gui.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        gui.setSize(1200,800);

        setupSim(); 

        JMenuBar bar = new JMenuBar();
        JMenu fileMenu = new JMenu("File");
        fileMenu.add(new JMenuItem("Open"));
        fileMenu.add(new JMenuItem("Save"));
        JMenu geometryMenu = new JMenu("Geometry");
        geometryMenu.add(new JMenuItem("Import"));
        geometryMenu.add(new JMenuItem("Add"));
        geometryMenu.add(new JMenuItem("Edit"));
        JMenu help = new JMenu("Help");
        bar.add(fileMenu);
        bar.add(geometryMenu);
        bar.add(help);

        gui.getContentPane().add(BorderLayout.NORTH, bar);

        gui.getContentPane().add(vp);

        gui.setVisible(true);

        sim.calc();

        vp.getResult = true;

        gui.repaint();
    }

    private static void plotTest() {
        int divisions = 100;
        Airfoil f = new Airfoil("23012",divisions);
        java.util.ArrayList<Double> x = new java.util.ArrayList<>();
        for(int i = 0; i <= divisions; i++) {
            x.add(i*1.0/divisions);
        }
        java.util.HashMap<String,Object> o = new java.util.HashMap<>();
        o.put("Title","Airfoil 23012");
        JPlot plot = new JPlot(x,f.getChamber(),o);

    }

    private static void setupSim() {
        Cartesian o = new Cartesian(1,0,0);

        Surface s1 = new Surface(o,0.5,1,5,5);
        Surface s2 = new Surface(o,0.5,-1,5,5);

        sim.surfaces.add(s1);
        sim.surfaces.add(s2);

        vp.setSurfaces(sim.surfaces);

        sim.setFreestream(50, 0.01, 0, 101325, 298);

        vp.setFreeStream(sim.getFreestream(),0.01);
        
    }

    private static void performanceTest() {
        long startTime , finishTime;
        int iter = 100000;

        startTime = System.nanoTime();
        for(int i = 0; i < iter; i++) {
            //double x = MyMath.quickCosd(2.1/iter);
            double x = MyMath.atan2(i/iter,i);
        }
        finishTime = System.nanoTime();
        System.out.println("Part 1 took: " + (finishTime - startTime)/1000 + " ms");

        startTime = System.nanoTime();
        for(int i = 0; i < iter; i++) {
            double x = Math.atan2(i/iter,i);
        }
        finishTime = System.nanoTime();
        System.out.println("Part 2 took: " + (finishTime - startTime)/1000 + " ms");

    }

    private static void valueTest() {
        SquareMatrix A = new SquareMatrix(5);
        A.set(0,0,2);
        A.set(0,1,-1);
        A.set(0,1,-1);
        A.set(0,1,-1);
        A.set(1,0,-1);
        A.set(1,1,2);
        A.set(1,2,-1);
        A.set(2,1,-1);
        A.set(2,2,2);
        A.set(2,3,-1);
        A.set(3,2,-1);
        A.set(3,3,2);

        Matrix b = new Matrix(5,1);
        b.set(0,0,1);
        b.set(1,0,8);
        b.set(2,0,-2);
        b.set(3,0,3);
        b.set(4,0,100);

        Matrix x = new Matrix(5,1);

        SquareMatrix.GaussSeidel(A,b,1e-6,x);

        System.out.println(x.get(0,0)+" "+x.get(1,0)+" "+x.get(2,0)+" "+x.get(3,0));

    }
    
}