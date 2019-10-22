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

    private Cartesian reference;

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

        AIC = new SquareMatrix(panels.size());
        b = new Matrix(panels.size(),1);
        // w = new Cartesian[panels.size()][panels.size()];

        for(int i = 0; i < panels.size();i++) {
            Panel p_i  = panels.get(i);
            Cartesian r_i = p_i.getCollocationPoint();
            Cartesian n_i = p_i.getNormal();
            for(int j = 0; j < panels.size(); j++) {
                AIC.set(i,j,panels.get(j).getInducedVelocityFactorAtPoint(r_i).dot(n_i));
            }
            b.set(i,0,freestream.dot(n_i));
        }

        
        Matrix x = new Matrix(panels.size(),1);

        SquareMatrix.SOR(AIC, b, 1e-12, x,0.01,0.2);

        // Matrix x = AIC.inv().mult(b);

        induced = AIC.mult(x);

        double LIFT = 0;
        double DRAG = 0;
        Cartesian MOMENT = new Cartesian();
        
        double rho = density*airspeed;
        for(int i = 0; i < panels.size(); i++) {
            Panel p_i = panels.get(i);
            // Cartesian v = new Cartesian(freestream);
            p_i.setCirculation(x.get(i,0));
            double circ = rho*x.get(i,0);

            double lift = circ*(p_i.vertices[2].y-p_i.vertices[0].y);

            double drag = Math.abs(circ*induced.get(i,0)*p_i.getNormal().z/airspeed);

            // Check this
            if(p_i.getNormal().z < 0) {
                LIFT += lift;
            } else {
                LIFT -= lift;
            }
            
            DRAG += drag;
            Cartesian d = freestream.mult(drag/airspeed);
            
            Cartesian l = p_i.getVortexVector().cross(freestream);
            l.normalize().multBy(lift);
            Cartesian f = new Cartesian(d.x+l.x,d.y+l.y,d.z+l.z);
            panels.get(i).setForce(f);
        }
        System.out.println("Lift = " + LIFT);
        System.out.println("Induced Drag = " + DRAG);
    }



}