package app;

import javax.swing.*;

import geometry.*;
import app.aero.VortexLatticeSteady;
import utils.Cartesian;
import utils.MyMath;

public class App {
    public static final String VERSION_NO = "0.2.0";

    static MainFrame gui;

    static VortexLatticeSteady sim = new VortexLatticeSteady();


    public static void main(String[] args) throws Exception {
        // valueTest();

        simTest();
    }

    private static void simTest() {
        gui = new MainFrame(VERSION_NO);

        gui.setSize(1200,800);
        gui.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        gui.setVisible(true);

        setupSim();
        sim.run();

        sim.printResults(1,0.5,"results_"+System.currentTimeMillis()+".dat");

        gui.repaint();
    }

    private static void setupSim() {
        /*
        sim.surfaces.add(new Surface(new Cartesian(1,0,0),0.5,1,8,8,5*MyMath.DEG2RAD,0,0));
        sim.surfaces.add(new Surface(new Cartesian(1,0,0),0.5,-1,8,8,-5*MyMath.DEG2RAD,0,0));
        */

        Airfoil NACA2402 = new Airfoil("2402",10);

        Cartesian root = new Cartesian(0,0,0);
        Cartesian tip = new Cartesian(0,1,0);
        Cartesian tip2 = new Cartesian(0,-1,0);
        sim.surfaces.add(new Surface(root,NACA2402,tip,NACA2402,0.5,0.5,-5*MyMath.DEG2RAD,0,6,10));
        sim.surfaces.add(new Surface(root,NACA2402,tip2,NACA2402,0.5,0.5,-5*MyMath.DEG2RAD,0,6,10));

        // Horizontal Tail
         Surface HStab = new Surface(new Cartesian(2.5,-0.5,0),0.25,1,4,6,0,0,0);
         HStab.addControlSurface(Surface.ELEVATOR, -0.4, 0.4, 0.3, 10);
         sim.surfaces.add(HStab);
        // sim.surfaces.add(new Surface(new Cartesian(2.5,-0.5,0),0.25,1,4,6,0,0,0));
        // Vertical Tail
        sim.surfaces.add(new Surface(new Cartesian(2.5,0,0),0.25,0.5,4,3,Math.PI/2,10*MyMath.DEG2RAD,5*MyMath.DEG2RAD));

        gui.vp.setSurfaces(sim.surfaces);

        sim.setFreestream(200, 0.01, 0, 101325, 298);

        sim.setReferencePoint(new Cartesian(0.5,0,0));

        gui.vp.setFreeStream(sim.getFreestream(),1.0/sim.getFreestream().r);
        
    }

    /*

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
        A.set(4,1,1);
        A.set(4,3,3);
        A.set(4,4,1);

        Matrix b = new Matrix(5,1);
        b.set(0,0,1);
        b.set(1,0,8);
        b.set(2,0,-2);
        b.set(3,0,3);
        b.set(4,0,100);

        Matrix x = new Matrix(5,1);

        SquareMatrix.GaussSeidel(A,b,1e-6,x);

        Matrix x2 = new Matrix(5,1);

        SquareMatrix.SOR(A,b,1e-9,x2,0.1,0.5);

        double det = A.det();

        SquareMatrix in = A.inv();

        Matrix x3 = new Matrix(5,1);

        try {
            Object[] arr = A.PLUDecompose();

            SquareMatrix.LUSolve((SquareMatrix)arr[1],(int[]) arr[0],b,x3);
        } catch (MatrixException e) {

        }

        System.out.println(x.get(0,0)+" "+x.get(1,0)+" "+x.get(2,0)+" "+x.get(3,0));

    }
    */
    

    
}