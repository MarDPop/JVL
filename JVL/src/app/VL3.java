package app;

import utils.*;

import geometry.*;

import java.util.ArrayList;
import java.util.Arrays;


/**
 * Vortex Lattice method, use http://www.cta-dlr2009.ita.br/Proceedings/PDF/59306.pdf for reference
 */
public class VL3 {

    public static double GAS_R = 8.31445; // J/ mol k

    public static double AIR_R = 287.1012; // J/ kg K

    private static double CONST_R = 401.9416; // J/ kg K

    private static double DYNAMIC_VISCOSITY_300 = 18.45e-6; //Pa s

    private static double fourpi = 0.25/Math.PI;

    /**
     * Free stream vector (m/s)
     */
    private Cartesian freestream = new Cartesian(); // X is considered nominal

    /**
     * Angle of Attack angle (rad)
     */
    private double alpha; 

    /**
     * Side slip angle (rad)
     */
    private double beta;

    /**
     * Static Air pressure (Pa)
     */
    private double staticPressure;

    /**
     * Airspeed (m/s) [ also magnitude of freestream vector]
     */
    private double airspeed; 

    /**
     * Freestream density (kg/m3)
     */
    private double density; 

    /**
     * Freestream temperature (K)
     */
    private double temperature;

    /**
     * Freestream Speed of sound (m/s)
     */
    private double a0;

    /**
     * Freestream mach
     */
    private double Mach;

    /**
     * Mach squared
     */
    private double M2;

    /**
     * Compressible Constant Beta 
     */
    private double MachConstant;

    /**
     * Freestream Reynolds Number factor (rho u/ mu) 1 / m
     */
    private double RE;

    /**
     * List of surfaces used in calculation
     */
    public ArrayList<Surface> surfaces = new ArrayList<Surface>();

    /**
     * List of panels contained in surfaces
     */
    private ArrayList<Panel> panels = new ArrayList<Panel>();

    /**
     * Array of vectors of induced velocity components of vortex based on geometry
     */
    private Cartesian[][] w;

    /**
     * Aerodynamic Influence Coefficient Matrix
     */
    private SquareMatrix AIC;

    /**
     * Vector to solve for (induced velocity normal to panel)
     */
    private Matrix b;

    /**
     * total force on vehicle
     */
    private Cartesian totalForce;

    /**
     * total moments on vehicle
     */
    private Cartesian totalMoment;

    /**
     * Reference point for moments
     */
    private Cartesian reference = new Cartesian();

    /**
     * Empty Constructor
     */
    public VL3() {}

    /**
     * Sets the surfaces
     * @param s
     */
    public void setSurfaces(ArrayList<Surface> s) {
        this.surfaces = s;
    }

    /**
     * Sets the freestream 
     * @param airspeed
     * @param alpha
     * @param beta
     * @param staticPressure
     * @param temperature
     */
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
        this.RE = density*airspeed/calcDynamicViscosity(temperature);
    }

    /**
     * Dynamic viscosity based on temperature curve fit
     * @param temperature
     * @return
     */
    private double calcDynamicViscosity(double temperature) {
        return temperature*(0.07549-4.643e-5*temperature);
    }

    /**
     * Gets dynamic pressure of freestream
     * @return
     */
    public double getQ() {
        return 0.7*this.staticPressure*this.M2;
    }

    /**
     * Calculation. 
     * Summary: Collects all the panels
     */
    public void calc() {
        // Collect Panels
        for(Surface s : surfaces) {
            for(Panel[] ps: s.getPanels()) {
                this.panels.addAll(Arrays.asList(ps));
            }
        }

        // Init Variables for Calculation
        int n = panels.size();
        AIC = new SquareMatrix(n);
        b = new Matrix(n,1);
        w = new Cartesian[panels.size()][panels.size()];

        for(int i = 0; i < n;i++) {
            // Panel being considered and some repeated vectors
            Panel p_i  = panels.get(i);
            Cartesian r_i = p_i.getCollocationPoint();
            Cartesian n_i = p_i.getNormal();
            for(int j = 0; j < n; j++) {
                // Get the contribution of panel j vortex on panel i
                w[i][j] = panels.get(j).getInducedVelocityFactorAtPoint(r_i);
                // get the only the normal component of contribution
                AIC.set(i,j,w[i][j].dot(n_i));
            }
            // set the sum of contributions to negate the freestream component normal to panel
            b.set(i,0,-freestream.dot(n_i));
        }

        // Solution Variable
        Matrix x = new Matrix(n,1);

        // Solve
        try {
            // Try LU factorization
            Object[] arr = AIC.PLUDecompose();
            SquareMatrix.LUSolve((SquareMatrix)arr[1],(int[]) arr[0],b,x);
            
        } catch (MatrixException e) {
            // Use Successive Overrelaxation if Fails
            System.out.println(e.getMessage());
            SquareMatrix.SOR(AIC, b, 1e-8, x,0.001,0.5);
        }

        // Set circulation solution to panels 
        for(int i = 0; i < panels.size(); i++) {
            panels.get(i).setCirculation(x.get(i,0));
        }

        // Init total forces/moments
        totalForce = new Cartesian();
        totalMoment = new Cartesian();

        // Compute forces and moments from panels
        for(int i = 0; i< n;i++){
            // Get panel considered
            Panel p_i = panels.get(i);
            Cartesian center = p_i.getCenter();
            // Vector for local induced velocity at center
            Cartesian local_induced = new Cartesian();
            for(int j = 0; j < n; j++){
                // panel j
                Panel p_j = panels.get(j);
                // Get induced velocity and multiply by circulation strength
                local_induced.addTo(p_j.getInducedVelocityFactorAtPoint(center).multBy(p_j.getCirculation()));
            }
            // Add to freestream velocity to get flow at center (would not necessarily be tangent)
            Cartesian local_velocity = freestream.add(local_induced);
            // Force computation is v x omega
            Cartesian local_force = local_velocity.cross(p_i.getVortexVector()).multBy(density*p_i.getCirculation());
            // Set force for panel
            p_i.setForce(local_force);
            // Sum across all panels
            totalForce.addTo(local_force);
            // Moment is moment arm x force
            totalMoment.addTo(local_force.cross(p_i.getCenter().sub(reference)));
        }

        // Compute Lift and Drag vectors from total force
        double LIFT = (-totalForce.x*freestream.z+totalForce.z*freestream.x)/airspeed;
        double DRAG = (totalForce.x*freestream.x+totalForce.z*freestream.z)/airspeed;

        System.out.println("Lift = " + LIFT + "N");
        System.out.println("Induced Drag = " + DRAG + "N");
        System.out.println("Pitching Moment = " + totalMoment.y + "N-m");
    }

    /* Getters and Setters */
    public Cartesian getFreestream() {
        return this.freestream;
    }

    public void setReferencePoint(Cartesian point) {
        this.reference = point;
    }

    public Cartesian getForces() {
        return totalForce;
    }

    public Cartesian getMoments() {
        return totalMoment;
    }

}