package app.aero;

import utils.*;

import geometry.*;

import java.util.ArrayList;
import java.util.Arrays;

public class VL {

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

    public VL() {}

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
        w = new Cartesian[panels.size()][panels.size()];

        for(int i = 0; i < panels.size();i++) {
            Panel p_i  = panels.get(i);
            Cartesian r_i = p_i.getCollocationPoint();
            double tmp;
            for(int j = 0; j < panels.size(); j++) {
                Panel p_j  = panels.get(j);

                Cartesian r1 = p_j.getA().sub(r_i);
                Cartesian r2 = p_j.getB().sub(r_i);

                Cartesian VAB = r1.cross(r2);
                VAB.normalize();
                double omega = p_j.getVortexVector().dot(r1)/r1.getMagnitude() - p_j.getVortexVector().dot(r2)/r2.getMagnitude();
                VAB.multBy(omega);

                Cartesian VA = new Cartesian();
                VA.y = -r1.z;
                VA.z = r1.y;
                tmp = r1.z*r1.z+r1.y*r1.y;
                tmp = (1 + r1.x/Math.sqrt(r1.x*r1.x + tmp))/tmp; 
                VA.multBy(tmp);

                Cartesian VB = new Cartesian();
                VB.y = r2.z;
                VB.z = -r2.y;
                tmp = r2.z*r2.z+r2.y*r2.y;
                tmp = (1 + r2.x/Math.sqrt(r2.x*r2.x + tmp))/tmp; 
                VB.multBy(tmp);

                w[i][j] = VAB.add(VA).add(VB);
                w[i][j].multBy(fourpi);

                AIC.set(i,j,w[i][j].dot(p_i.getNormal()));
            }
            b.set(i,0,-freestream.dot(p_i.getNormal()));
        }

        Matrix x = new Matrix(panels.size(),1);

        SquareMatrix.GaussSeidel(AIC, b, 1e-6, x);

        induced = AIC.mult(x);

        double LIFT = 0;
        double DRAG = 0;
        Cartesian MOMENT = new Cartesian();
        
        for(int i = 0; i < panels.size(); i++) {
            Panel p_i = panels.get(i);
            // Cartesian v = new Cartesian(freestream);
            p_i.setCirculation(b.get(i,0));
            double rho = density*b.get(i,0);
            double lift = rho*airspeed*Math.abs(p_i.vertices[2].y-p_i.vertices[0].y);
            LIFT += lift;
            double drag = rho*induced.get(i,0)/airspeed;
            Cartesian d = freestream.mult(drag);
            DRAG += drag;
            Cartesian l = p_i.getVortexVector().cross(freestream);
            l.normalize().multBy(lift);
            Cartesian f = new Cartesian(d.x+l.x,d.y+l.x,d.z+l.x);
            panels.get(i).setForce(f);
        }
        System.out.println("Lift = " + LIFT);
        System.out.println("Induced Drag = " + DRAG);
    }



}