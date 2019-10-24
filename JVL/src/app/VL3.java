package app;

import utils.*;

import geometry.*;

import java.util.ArrayList;
import java.util.Arrays;


/**
 * Vortex Lattice method, use http://www.cta-dlr2009.ita.br/Proceedings/PDF/59306.pdf for reference
 */
public class VL3 {

    public static double GAS_R = 8.31445;

    public static double AIR_R = 287.1012;

    private static double CONST_R = 401.9416;


    private Cartesian freestream = new Cartesian(); // X is considered nominal

    private double alpha; // angle of attack

    private double beta; // sideslip

    private double staticPressure;

    private double airspeed; // airspeed

    private double density; // air density (not needed for calc)

    private double temperature;

    private double a0;

    private double Mach;

    private double M2;

    private double MachConstant;

    private double fourpi = 0.25/Math.PI;

    public ArrayList<Surface> surfaces = new ArrayList<Surface>();

    private ArrayList<Panel> panels = new ArrayList<Panel>();

    private Cartesian[][] w;

    private Matrix induced;

    private Cartesian[] Force;

    private SquareMatrix AIC;

    private Matrix b;

    private Cartesian totalForce;

    private Cartesian totalMoment;

    private Cartesian reference = new Cartesian();

    public VL3() {}

    public void setSurfaces(ArrayList<Surface> s) {
        this.surfaces = s;
    }

    public void setFreestream(double airspeed, double alpha, double beta, double staticPressure, double temperature) {
        freestream.x = airspeed*Math.cos(alpha)*Math.cos(beta);
        freestream.y = -airspeed*Math.sin(beta);
        freestream.z = airspeed*Math.sin(alpha)*Math.cos(beta);
        freestream.r = airspeed;

        this.airspeed = airspeed;
        this.alpha = alpha;
        this.beta = beta;

        this.staticPressure = staticPressure;
        this.density = staticPressure/(AIR_R*temperature);
        this.temperature = temperature;

        this.a0 = Math.sqrt(CONST_R*temperature);

        this.Mach = airspeed/a0;
        this.M2 = Mach*Mach;
        this.MachConstant = 1-M2;
    }

    public double getQ() {
        return 0.7*this.staticPressure*this.M2;
    }

    public Cartesian getFreestream() {
        return this.freestream;
    }

    public void setReferencePoint(Cartesian point) {
        this.reference = point;
    }

    public void calc() {
        for(Surface s : surfaces) {
            for(Panel[] ps: s.getPanels()) {
                this.panels.addAll(Arrays.asList(ps));
            }
        }

        int n = panels.size();
        AIC = new SquareMatrix(n);
        b = new Matrix(n,1);
        w = new Cartesian[panels.size()][panels.size()];

        for(int i = 0; i < n;i++) {
            Panel p_i  = panels.get(i);
            Cartesian r_i = p_i.getCollocationPoint();
            Cartesian n_i = p_i.getNormal();
            for(int j = 0; j < n; j++) {
                w[i][j] = panels.get(j).getInducedVelocityFactorAtPoint(r_i);
                AIC.set(i,j,w[i][j].dot(n_i));
            }
            b.set(i,0,-freestream.dot(n_i));
        }

        Matrix x = new Matrix(n,1);

        try {
            Object[] arr = AIC.PLUDecompose();
            SquareMatrix.LUSolve((SquareMatrix)arr[1],(int[]) arr[0],b,x);
            
        } catch (MatrixException e) {
            System.out.println(e.getMessage());
            SquareMatrix.SOR(AIC, b, 1e-8, x,0.001,0.5);
        }

        induced = AIC.mult(x);

        for(int i = 0; i < panels.size(); i++) {
            panels.get(i).setCirculation(x.get(i,0));
        }

        totalForce = new Cartesian();
        totalMoment = new Cartesian();

        for(int i = 0; i< n;i++){
            Panel p_i = panels.get(i);
            Cartesian center = p_i.getCenter();
            Cartesian local_induced = new Cartesian();
            for(int j = 0; j < n; j++){
                Panel p_j = panels.get(j);
                local_induced.addTo(p_j.getInducedVelocityFactorAtPoint(center).multBy(p_j.getCirculation()));
            }
            Cartesian local_velocity = freestream.add(local_induced);
            Cartesian local_force = local_velocity.cross(p_i.getVortexVector()).multBy(density*p_i.getCirculation());
            p_i.setForce(local_force);
            totalForce.addTo(local_force);
            totalMoment.addTo(local_force.cross(p_i.getCenter().sub(reference)));
        }

        double LIFT = (-totalForce.x*freestream.z+totalForce.z*freestream.x)/airspeed;
        double DRAG = (totalForce.x*freestream.x+totalForce.z*freestream.z)/airspeed;

        System.out.println("Lift = " + LIFT);
        System.out.println("Induced Drag = " + DRAG);
    }



}