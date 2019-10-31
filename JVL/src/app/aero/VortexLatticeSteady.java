package app.aero;

import utils.*;
import app.Options;
import geometry.*;

import java.util.ArrayList;
import java.util.Arrays;

import files.ImportExport;


/**
 * Vortex Lattice method, use http://www.cta-dlr2009.ita.br/Proceedings/PDF/59306.pdf for reference
 */
public class VortexLatticeSteady {

    public static double GAS_R = 8.31445; // J/ mol k

    public static double AIR_R = 287.1012; // J/ kg K

    private static double CONST_R = 401.9416; // J/ kg K

    // private static double DYNAMIC_VISCOSITY_300 = 18.45e-6; //Pa s

    public Options options = new Options();

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
    private double staticTemperature;

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
     * Potential Flow Scaling factors
     */
    private Cartesian compressibleScaling;

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
     * 0 = dForce/dAlpha 1 = dMoment/dAlpha 2= dForce/dBeta 3 = dMoment/dBeta 4 = dForce/du 5 = dMoment/du 6 = dForce/dAirspeed 7 = dMoment/dAirspeed
     */
    private Cartesian[] derivatives = new Cartesian[8];

    /**
     * Vehicle lift
     */
    private double VehicleLift;

    /**
     * Vehicle Induced Drag
     */
    private double VehicleInducedDrag;

    /**
     * Vehicle Drag
     */
    private double VehicleDrag;

    /**
     * number of panels
     */
    int n;

    /**
     * Empty Constructor
     */
    public VortexLatticeSteady() {
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
        this.staticTemperature = temperature;
        this.density = staticPressure/(AIR_R*this.staticTemperature);

        // Speed of sound
        this.a0 = Math.sqrt(CONST_R*temperature);

        // Mach Number and Constants associated
        this.Mach = airspeed/a0;
        this.M2 = Mach*Mach;
        this.MachConstant = Math.sqrt(1-M2);

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
     * Gets compressible density ratio
     * @return
     */
    public double getCompressibleDensityRatio(double Mach) {
        double Beta = 1-0.8*Mach*Mach;
        return Beta*Beta*Math.sqrt(Beta);
    }

    /**
     * Calculation. 
     */
    public void calc(boolean compressibility) {
        if ( n == 0 )
            initCalc();

        if (compressibility) {
            calcCompressible();
        } else {
            calcIncompressible();
        }
        getForcesAndMoments();
    }

    /**
     * Initializes variables for calculation
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
     * 
     */
    private void calcCompressible() {
        calcIncompressible();
        compressibleScaling = new Cartesian();
        
        // Solve
        try {
            // Try LU factorization
            int count = 0;
            double err = 1;
            while (err > 1e-3 && count < 10){
                computeAICCompressible();
                Object[] arr = AIC.PLUDecompose();
                SquareMatrix.LUSolve((SquareMatrix)arr[1],(int[]) arr[0],b,x);
                
            }
        } catch (MatrixException e) {
            // Use Successive Overrelaxation if Fails
            System.out.println(e.getMessage());
            computeAICCompressible();
            SquareMatrix.SOR(AIC, b, 1e-8, x,0.001,0.5);
        }
    }

    private void computeAICCompressible() {
        // Define AIC and velocity to cancel out
        for(int i = 0; i < n;i++) {
            // Panel being considered and some useful repeated vectors
            Panel p_i  = panels.get(i);
            Cartesian r_i = p_i.getCollocationPoint();
            Cartesian n_i = p_i.getCompressibleNormal();
            for(int j = 0; j < n; j++) {
                // Get the contribution of panel j vortex on panel i
                w[i][j] = panels.get(j).getInducedVelocityFactorAtPointCompressible(r_i,compressibleScaling);
                // get the only the normal component of contribution
                AIC.set(i,j,w[i][j].dot(n_i));
            }
            // set the sum of contributions to negate the freestream component normal to panel
            b.set(i,0,-freestream.dot(n_i));
        }
        // Set circulation solution to panels 
        for(int i = 0; i < panels.size(); i++) {
            panels.get(i).setCirculation(x.get(i,0));
        }
    }

     /**
     * 
     */
    private void calcIncompressible() {
        // Define AIC and velocity to cancel out
        for(int i = 0; i < n;i++) {
            // Panel being considered and some useful repeated vectors
            Panel p_i  = panels.get(i);
            Cartesian r_i = p_i.getCollocationPoint();
            Cartesian n_i = p_i.getNormal();
            for(int j = 0; j < n; j++) {
                // Get the contribution of panel j vortex on panel i
                w[i][j] = panels.get(j).getInducedVelocityFactorAtPointIncompressible(r_i);
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
        } catch (MatrixException e) {
            // Use Successive Overrelaxation if Fails
            System.out.println(e.getMessage());
            SquareMatrix.SOR(AIC, b, 1e-8, x,0.001,0.5);
        }

        // Set circulation solution to panels 
        for(int i = 0; i < panels.size(); i++) {
            panels.get(i).setCirculation(x.get(i,0));
        }
    }

    /**
     * 
     */
    private void getForcesAndMoments() {

        // Init total forces/moments
        totalForce = new Cartesian();
        totalMoment = new Cartesian();

        double scale = 1/MachConstant;
        
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
                if (options.compressible) {

                } else {
                    local_induced.addTo(p_j.getInducedVelocityFactorAtPointIncompressible(center).multBy(p_j.getCirculation()));
                }
                
            }
            // Add to freestream velocity to get flow at center (would not necessarily be tangent)
            Cartesian local_velocity = freestream.add(local_induced);
            // Force computation is v x omega
            Cartesian local_force = local_velocity.cross(p_i.getVortexVector()).multBy(density*p_i.getCirculation());
            if(!options.compressible && options.scaleCompressible) {
                local_force.multBy(scale);
            }
            // Set force for panel
            p_i.setForce(local_force);
            // Sum across all panels
            totalForce.addTo(local_force);
            // Moment is moment arm x force
            totalMoment.addTo(p_i.getCenter().sub(reference).cross(local_force));
        }

        // Compute Lift and Drag vectors from total force
        VehicleLift = (-totalForce.x*freestream.z+totalForce.z*freestream.x)/airspeed;
        VehicleInducedDrag = (totalForce.x*freestream.x+totalForce.y*freestream.y+totalForce.z*freestream.z)/airspeed;
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
        VehicleDrag = VehicleInducedDrag + drag;
        return drag;
    }

    /**
     * Runs the Vortex Lattice for nominal condition and deviated positions to get approximate stability derivatives
     * @param refArea
     * @param refLength
     */
    public void run() {
        // Calculate initial condition
        calc(options.compressible);
        // Save Forces and Momoments
        Cartesian NomForces = new Cartesian(totalForce);
        Cartesian NomMoments = new Cartesian(totalMoment);
        double L = VehicleLift;
        double D = VehicleInducedDrag;

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

        calc(options.compressible);
        F = new Cartesian(totalForce);
        M = new Cartesian(totalMoment);

        freestream.x = airspeed*Math.cos(alpha+dAlpha)*Math.cos(beta);
        freestream.y = airspeed*Math.sin(beta);
        freestream.z = airspeed*Math.sin(alpha+dAlpha)*Math.cos(beta);
        calc(options.compressible);

        this.derivatives[0] = totalForce.sub(F).mult(0.5/dAlpha);
        this.derivatives[1]  = totalMoment.sub(M).mult(0.5/dAlpha);

        // Beta
        freestream.x = airspeed*Math.cos(dAlpha)*Math.cos(beta-dBeta);
        freestream.y = airspeed*Math.sin(beta-dBeta);
        freestream.z = airspeed*Math.sin(dAlpha)*Math.cos(beta-dBeta);

        calc(options.compressible);
        F = new Cartesian(totalForce);
        M = new Cartesian(totalMoment);

        freestream.x = airspeed*Math.cos(alpha+dAlpha)*Math.cos(beta+dBeta);
        freestream.y = airspeed*Math.sin(beta+dBeta);
        freestream.z = airspeed*Math.sin(alpha+dAlpha)*Math.cos(beta+dBeta);
        calc(options.compressible);
        
        this.derivatives[2] = totalForce.sub(F).mult(0.5/dBeta);
        this.derivatives[3]  = totalMoment.sub(M).mult(0.5/dBeta);

        // Acceleration in X axis
        freestream.x = (airspeed-dAirspeed)*Math.cos(alpha)*Math.cos(beta);
        freestream.y = airspeed*Math.sin(beta);
        freestream.z = airspeed*Math.sin(alpha)*Math.cos(beta);

        calc(options.compressible);

        F = new Cartesian(totalForce);
        M = new Cartesian(totalMoment);

        freestream.x = (airspeed+dAirspeed)*Math.cos(alpha)*Math.cos(beta);
        freestream.y = airspeed*Math.sin(beta);
        freestream.z = airspeed*Math.sin(alpha)*Math.cos(beta);
        calc(options.compressible);
        
        this.derivatives[4] = totalForce.sub(F).mult(0.5/dAirspeed);
        this.derivatives[5]  = totalMoment.sub(M).mult(0.5/dAirspeed);

        // Airspeed
        freestream.x = (airspeed-dAirspeed)*Math.cos(alpha)*Math.cos(beta);
        freestream.y = (airspeed-dAirspeed)*Math.sin(beta);
        freestream.z = (airspeed-dAirspeed)*Math.sin(alpha)*Math.cos(beta);

        calc(options.compressible);

        F = new Cartesian(totalForce);
        M = new Cartesian(totalMoment);

        freestream.x = (airspeed+dAirspeed)*Math.cos(alpha)*Math.cos(beta);
        freestream.y = (airspeed+dAirspeed)*Math.sin(beta);
        freestream.z = (airspeed+dAirspeed)*Math.sin(alpha)*Math.cos(beta);
        calc(options.compressible);
        
        this.derivatives[6] = totalForce.sub(F).mult(0.5/dAirspeed);
        this.derivatives[7]  = totalMoment.sub(M).mult(0.5/dAirspeed);

        // Note: typically speed derivatives are in component direction but airspeed derivative has merit for wind change
        
        // Reset to base
        this.totalForce = NomForces;
        this.totalMoment = NomMoments;
        this.VehicleLift = L;
        this.VehicleInducedDrag = D;

        // Approximate skin drag
        calcApproxSkinDrag(); 
        // NOTE: just assume that skin drag scales with velocity squared for dCD/du (ie 2*u*drag_skin/v)
    }

    /**
     * Prints results 
     */
    public void printResults(double refArea, double refLength, String filename) {
        double refForce = this.getQ()*refArea;
        double refMoment = refForce*refLength;

        // Print a few to screen
        System.out.println("Lift = " + VehicleLift + "N");
        System.out.println("Induced Drag = " + VehicleInducedDrag + "N");
        System.out.println("Total Drag = " + VehicleDrag +"N");
        System.out.println("Pitching Moment = " + totalMoment.y +"N-m");
        System.out.println("Cm_alpha = " + derivatives[1].y +"N-m/rad");

        // Print all to file
        // Collect Coeffients
        String[] coefNames= new String[7];
        double[] coefficients = new double[7];
        coefNames[0] = "C_L (Lift)";
        coefficients[0] = VehicleLift/refForce;
        coefNames[1] = "C_D (Drag)";
        coefficients[1] = VehicleDrag/refForce;
        coefNames[2] = "C_M (Pitching Moment)";
        coefficients[2] = totalMoment.y/refMoment;
        coefNames[3] = "C_M_alpha (Pitch damping)";
        coefficients[3] = derivatives[1].y/refMoment;
        coefNames[4] = "C_N (Yaw Moment)";
        coefficients[4] = totalMoment.z/refMoment;
        coefNames[5] = "C_L (Roll Moment)";
        coefficients[5] = totalMoment.x/refMoment;
        coefNames[6] = "C_N_beta (Roll damping)";
        coefficients[6] = derivatives[3].z/refMoment;

        // Other variables
        double[] references = new double[]{airspeed, alpha, beta, Mach, staticPressure, refArea, refLength, reference.x,reference.y,reference.z};
        
        // Print
        ImportExport.writeCoefficients(coefficients, coefNames, references, filename);
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

    public double getLift() {
        return VehicleLift;
    }

}