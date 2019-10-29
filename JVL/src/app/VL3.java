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
     * Vector of solution (Circulation Strengths)
     */
    private Matrix x;

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
     *  All coefficients (force, moments, stability) 
     */
    private double[][] coefficents = new double[7][3]; // CFx, CFy, CFz; CMx, CMy, CMz; CMy_alpha, CMz_beta

    /**
     * 0 = dForce/dAlpha 1 = dMoment/dAlpha 2= dForce/dBeta 3 = dMoment/dBeta 4 = dForce/dAirspeed 5 = dMoment/dAirspeed
     */
    private Cartesian[] derivatives = new Cartesian[6];

    /**
     * Vehicle lift
     */
    private double LIFT;

    /**
     * Vehicle Drag
     */
    private double DRAG;

    /**
     * number of panels
     */
    int n;

    /**
     * Empty Constructor
     */
    public VL3() {
        n = 0;
    }

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
        // setup free stream vector from angle of attack and side slip
        freestream.x = airspeed*Math.cos(alpha)*Math.cos(beta);
        freestream.y = airspeed*Math.sin(beta);
        freestream.z = airspeed*Math.sin(alpha)*Math.cos(beta);
        freestream.r = airspeed;

        // Definition of Free Stream Vector stored
        this.airspeed = airspeed;
        this.alpha = alpha;
        this.beta = beta;

        // Freestream Air State
        this.staticPressure = staticPressure;
        this.temperature = temperature;
        this.density = staticPressure/(AIR_R*temperature);

        // Speed of sound
        this.a0 = Math.sqrt(CONST_R*temperature);

        // Mach Number and Constants associated
        this.Mach = airspeed/a0;
        this.M2 = Mach*Mach;
        this.MachConstant = 1-M2;

        // Reynolds number factor
        this.RE = density*airspeed/calcDynamicViscosity(temperature);
    }

    /**
     * Dynamic viscosity based on temperature curve fit
     * @param temperature
     * @return dynamic viscosity (Pa s)
     */
    private double calcDynamicViscosity(double temperature) {
        return temperature*(0.07549-4.643e-5*temperature)*1e-6;
    }

    /**
     * Gets dynamic pressure of freestream
     * @return dynamic pressure in Pa
     */
    public double getQ() {
        return 0.7*this.staticPressure*this.M2;
    }

    /**
     * 
     */
    private void initCalc() {
        // Collect Panels
        for(Surface s : surfaces) {
            for(Panel[] ps: s.getPanels()) {
                this.panels.addAll(Arrays.asList(ps));
            }
        }

        // Init Variables for Calculation
        n = panels.size();
        // Aerodynamic Influence Coefficent Matrix
        AIC = new SquareMatrix(n);
        // Solution Variable (flow to cancel out)
        b = new Matrix(n,1);
        // Induced velocity vector
        w = new Cartesian[panels.size()][panels.size()];
        // Solution Variable
        x = new Matrix(n,1);
    }

    /**
     * Calculation. 
     */
    public void calc() {
        if ( n == 0 )
            initCalc();

        for(int i = 0; i < n;i++) {
            // Panel being considered and some useful repeated vectors
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

        // Solve
        try {
            // Try LU factorization
            Object[] arr = AIC.PLUDecompose();
            SquareMatrix.LUSolve((SquareMatrix)arr[1],(int[]) arr[0],b,x);
            System.out.println("solution found LU");
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
        LIFT = (-totalForce.x*freestream.z+totalForce.z*freestream.x)/airspeed;
        DRAG = (totalForce.x*freestream.x+totalForce.z*freestream.z)/airspeed;
    }

    /**
     * Calculates approximate skin drag from panels and a reynolds number function
     * @return
     */
    public double calcApproxSkinDrag() {
        // Skin Friction coefficient, reynolds number
        double Cf, reynolds;
        // Dynamic Pressure (use a factor of 2 for top and bottom surface)
        double q = 2*this.getQ();
        // Drag
        double drag = 0;
        for(Surface s : surfaces) {
            for(Panel[] ps: s.getPanels()) {
                // Leading edge x used to calculate distance for reynolds number
                double LE_x = ps[0].vertices[0].x;
                for (Panel p : ps) {
                    // reynolds number approx
                    reynolds = this.RE*(p.getCenter().x-LE_x);
                    if (reynolds < 3000) {
                        // Blasius
                        Cf = 0.664/Math.sqrt(reynolds);
                    } else {
                        // Prandtl
                        Cf = 0.027/Math.pow(reynolds,0.145);
                    }
                    // sum drag
                    drag += Cf*q*p.getArea();
                }
            }
        }
        return drag;
    }

    public void run(double refArea, double refLength) {
        // Calculate initial condition
        calc();
        // Save Forces and Momoments
        Cartesian NomForces = new Cartesian(totalForce);
        Cartesian NomMoments = new Cartesian(totalMoment);
        double L = LIFT;
        double D = DRAG;

        // Temp Force and Moment variables
        Cartesian F, M;

        // Discrete velocity changes
        double dAlpha = 0.01;
        double dBeta = 0.01;
        double dAirspeed = 1;

        // Alpha
        freestream.x = airspeed*Math.cos(alpha-dAlpha)*Math.cos(beta);
        freestream.y = airspeed*Math.sin(beta);
        freestream.z = airspeed*Math.sin(alpha-dAlpha)*Math.cos(beta);

        calc();
        F = new Cartesian(totalForce);
        M = new Cartesian(totalMoment);

        freestream.x = airspeed*Math.cos(alpha+dAlpha)*Math.cos(beta);
        freestream.y = airspeed*Math.sin(beta);
        freestream.z = airspeed*Math.sin(alpha+dAlpha)*Math.cos(beta);
        calc();

        this.derivatives[0] = totalForce.sub(F).mult(0.5/dAlpha);
        this.derivatives[1]  = totalMoment.sub(M).mult(0.5/dAlpha);

        // Beta
        freestream.x = airspeed*Math.cos(dAlpha)*Math.cos(beta-dBeta);
        freestream.y = airspeed*Math.sin(beta-dBeta);
        freestream.z = airspeed*Math.sin(dAlpha)*Math.cos(beta-dBeta);

        calc();
        F = new Cartesian(totalForce);
        M = new Cartesian(totalMoment);

        freestream.x = airspeed*Math.cos(alpha+dAlpha)*Math.cos(beta+dBeta);
        freestream.y = airspeed*Math.sin(beta+dBeta);
        freestream.z = airspeed*Math.sin(alpha+dAlpha)*Math.cos(beta+dBeta);
        calc();
        
        this.derivatives[2] = totalForce.sub(F).mult(0.5/dBeta);
        this.derivatives[3]  = totalMoment.sub(M).mult(0.5/dBeta);

        // Airspeed
        freestream.x = (airspeed-dAirspeed)*Math.cos(alpha)*Math.cos(beta);
        freestream.y = (airspeed-dAirspeed)*Math.sin(beta);
        freestream.z = (airspeed-dAirspeed)*Math.sin(alpha)*Math.cos(beta);

        calc();
        F = new Cartesian(totalForce);
        M = new Cartesian(totalMoment);

        freestream.x = (airspeed+dAirspeed)*Math.cos(alpha)*Math.cos(beta);
        freestream.y = (airspeed+dAirspeed)*Math.sin(beta);
        freestream.z = (airspeed+dAirspeed)*Math.sin(alpha)*Math.cos(beta);
        calc();
        
        this.derivatives[4] = totalForce.sub(F).mult(0.5/dBeta);
        this.derivatives[5]  = totalMoment.sub(M).mult(0.5/dBeta);

        // Constant for Coefficients
        double constant4Force = this.getQ()*refArea;
        double constant4Moment = constant4Force*refLength;
        
        // Reset to base
        this.totalForce = NomForces;
        this.totalMoment = NomMoments;
        this.LIFT = L;
        this.DRAG = D;

        // Approximate skin drag
        double drag = calcApproxSkinDrag();

        // Print out
        System.out.println("Lift = " + LIFT + "N");
        System.out.println("Induced Drag = " + DRAG + "N");
        System.out.println("Skin Drag = " + drag +"N");
        System.out.println("Pitching Moment = " + NomMoments.y +"N-m");
        System.out.println("Cm_alpha = " + derivatives[1].y +"N-m/rad");
    }



    /* Getters and Setters */
    public Cartesian getFreestream() {
        return this.freestream;
    }

    public void setReferencePoint(Cartesian point) {
        this.reference = point;
    }

    public Cartesian getReferencePoint() {
        return this.reference;
    }

    public Cartesian getForces() {
        return totalForce;
    }

    public Cartesian getMoments() {
        return totalMoment;
    }

}