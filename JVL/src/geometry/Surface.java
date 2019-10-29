package geometry;

import utils.Cartesian;
import utils.CartesianMatrix;
import utils.MyMath;

import java.util.ArrayList;

public class Surface {

    String name;
    
    Panel[][] panels;

    int nChord;

    int nSpan;

    double span;

    double root_chord;

    double taperRatio;

    double sweep;

    double dihedral;

    Cartesian origin = new Cartesian();

    Cartesian lowerBounds = new Cartesian();

    Cartesian upperBounds = new Cartesian();

    public static final int AILERON = 1;

    public static final int ELEVATOR = 2;

    public static final int RUDDER = 3;

    public static final int FLAP = 4;

    public static final int SPOILER = 5;

    // Control Surface array = 0: type 1: iStart 2:iEnd 3:jStart 4:jEnd 5: rotation in deg
    private ArrayList<int[]> controlSurface = new ArrayList<>();

    public Surface() {}

    /**
     * Load flat plate
     * @param origin
     * @param root_chord
     * @param span
     * @param nChord
     * @param nSpan
     * @param angle
     * @param taper
     * @param sweep
     */
    public Surface(Cartesian origin, double root_chord, double span,int nChord, int nSpan, double angle, double taperAngle, double sweep) {

        if(nChord < 0 || nSpan < 0 ) {
            nChord = Math.abs(nChord);
            nSpan = Math.abs(nSpan);
            System.out.println("Input number of divisions must be positive, absolute value used.");
        }

        if(root_chord < 0) {
            root_chord = Math.abs(root_chord);
            origin.x += root_chord;
            System.out.println("Chord must be positive, shifted origin upstream.");
        }

        this.nChord = nChord;
        this.nSpan = nSpan;
        this.span = span;
        this.root_chord = root_chord;
        this.taperRatio = 1-2*Math.tan(taperAngle)*span/root_chord;
        this.sweep = sweep;
        this.dihedral = angle;

        panels = new Panel[nSpan][nChord];
        
        this.origin = origin;
        this.lowerBounds.x = -1e+200;
        this.upperBounds.x = 1e+200;
        this.lowerBounds.y = -1e+200;
        this.upperBounds.y = 1e+200;
        this.lowerBounds.z = -1e+200;
        this.upperBounds.z = 1e+200;

        double d_span = span/nSpan;
        double deltaXSweep = Math.tan(sweep)*d_span;
        double spanRatioTaper = Math.tan(taperAngle)*d_span;
        double deltaXTaper = spanRatioTaper*0.5;

        double dy = d_span*Math.cos(angle);
        double dz = d_span*Math.sin(angle);

        for(int i = 0; i < nSpan; i++) {
            int i_outboard = i+1;
            double start_x_inboard = origin.x + i*(deltaXTaper + deltaXSweep);
            double chord_inboard = root_chord - i*spanRatioTaper;
            double dx_inboard = chord_inboard/nChord;
            double start_x_outboard = origin.x + i_outboard*(deltaXTaper + deltaXSweep);
            double chord_outboard = root_chord - i_outboard*spanRatioTaper;
            double dx_outboard = chord_outboard/nChord;
            for(int j = 0; j < nChord; j++) {
                panels[i][j] = new Panel();
                Cartesian x1y1 = new Cartesian(start_x_inboard + j*dx_inboard,origin.y + i*dy,i*dz);
                Cartesian x2y1 = new Cartesian(start_x_inboard +(j+1)*dx_inboard,origin.y +i*dy,i*dz);
                Cartesian x2y2 = new Cartesian(start_x_outboard +(j+1)*dx_outboard,origin.y +i_outboard*dy,i_outboard*dz);
                Cartesian x1y2 = new Cartesian(start_x_outboard +j*dx_outboard,origin.y +i_outboard*dy,i_outboard*dz);
                panels[i][j].setVertices(x1y1, x2y1, x2y2, x1y2);
            }
        }

    
    }

    /**
     * Load Surface from Tip and Root airfoils
     * @param LE1
     * @param root
     * @param LE2
     * @param tip
     * @param root_chord
     * @param tip_chord
     * @param root_twist
     * @param tip_twist
     * @param nSpan
     * @param nChord
     */
    public Surface(Cartesian LE1, Airfoil root, Cartesian LE2, Airfoil tip, double root_chord, double tip_chord, double root_twist, double tip_twist, int nSpan, int nChord) {
        
        if(root.chamber.size() == tip.chamber.size()) {
            int n = root.chamber.size();
            if(n == nChord+1) {
                // init panels
                this.nChord = root.chamber.size();
                this.nSpan = nSpan;

                panels = new Panel[nSpan][nChord];

                // Get leading edge
                Cartesian LE = LE2.sub(LE1);
                double dy = LE.y/nSpan;

                // Get definition
                this.span = LE.getMagnitude();
                this.root_chord = root_chord;
                this.taperRatio = tip_chord/root_chord;
                this.sweep = Math.atan((LE2.x+tip_chord*0.5-(LE1.x+root_chord*0.5))/span) ;
                this.dihedral = Math.asin(LE.z/span);

                // rotate and scale airfoil
                
                double[] x_root = new double[n];
                double[] z_root = new double[n];
                double[] x_tip = new double[n];
                double[] z_tip = new double[n];
                double c_root = Math.cos(root_twist);
                double s_root = Math.sin(root_twist);
                double c_tip = Math.cos(tip_twist);
                double s_tip = Math.sin(tip_twist);
                for (int i = 0; i < n; i++) {
                    double x = root_chord*i/nChord;
                    z_root[i] = root_chord*root.chamber.get(i);
                    x_root[i] = c_root*x - s_root*z_root[i];
                    z_root[i] = s_root*x + c_root*z_root[i];

                    x = tip_chord*i/nChord;
                    z_tip[i] = tip_chord*tip.chamber.get(i);
                    x_tip[i] = c_tip*x - s_tip*z_tip[i];
                    z_tip[i] = s_tip*x + c_tip*z_tip[i];
                }

                for(int i = 0; i < nSpan; i++) {
                    int i_outboard = i+1;
                    
                    double le_y_inboard = LE1.y + i*dy;
                    double le_y_outboard = LE1.y + i_outboard*dy;

                    for(int j = 0; j < nChord; j++) {
                        int j_te = j+1;
                        panels[i][j] = new Panel();
                        double dz_le = (z_tip[j]-z_root[j])/nSpan;
                        double dz_te = (z_tip[j_te]-z_root[j_te])/nSpan;
                        double dx_le = (x_tip[j]-x_root[j])/nSpan;
                        double dx_te = (x_tip[j_te]-x_root[j_te])/nSpan;
                        Cartesian x1y1 = new Cartesian(x_root[j] + i*dx_le, le_y_inboard, z_root[j]+i*dz_le);
                        Cartesian x2y1 = new Cartesian(x_root[j_te] + i*dx_te, le_y_inboard, z_root[j_te] + i*dz_te);
                        Cartesian x2y2 = new Cartesian(x_root[j_te] + i_outboard*dx_te, le_y_outboard, z_root[j_te] + i_outboard*dz_te);
                        Cartesian x1y2 = new Cartesian(x_root[j] +i_outboard*dx_le, le_y_outboard, z_root[j] + i_outboard*dz_le);
                        panels[i][j].setVertices(x1y1, x2y1, x2y2, x1y2);
                    }
                }
            } else {
                
            }
        }
    }

    /**
     * Load Surface from Tip and Root airfoils
     * @param LE1
     * @param root
     * @param LE2
     * @param tip
     * @param root_chord
     * @param tip_chord
     * @param root_twist
     * @param tip_twist
     * @param nSpan
     * @param nChord
     */
    public Surface(Airfoil root, Airfoil tip, double root_chord, double tip_chord, Cartesian origin, double span, double dihedral, double LEangle, double root_twist, double tip_twist, int nSpan, int nChord) {
        
        if(root.chamber.size() == tip.chamber.size()) {
            int n = root.chamber.size();
            if(n == nChord+1) {
                // init panels
                this.nChord = root.chamber.size();
                this.nSpan = nSpan;

                panels = new Panel[nSpan][nChord];

                // Get leading edge
                double dy = span/nSpan;

                // Get definition
                this.span = span;
                this.root_chord = root_chord;
                double x_LE = span*Math.sin(LEangle);
                this.taperRatio = tip_chord/root_chord;
                this.sweep = Math.asin((tip_chord*0.5 + x_LE - root_chord*0.5)/span); // sweep from mid chord to mid chord
                this.dihedral = dihedral;

                // rotate and scale airfoil
                
                double[] x_root = new double[n];
                double[] z_root = new double[n];
                double[] x_tip = new double[n];
                double[] z_tip = new double[n];
                double c_root = Math.cos(root_twist);
                double s_root = Math.sin(root_twist);
                double c_tip = Math.cos(tip_twist);
                double s_tip = Math.sin(tip_twist);
                for (int i = 0; i < n; i++) {
                    double x = root_chord*i/nChord;
                    z_root[i] = root_chord*root.chamber.get(i);
                    x_root[i] = c_root*x - s_root*z_root[i];
                    z_root[i] = s_root*x + c_root*z_root[i];

                    x = tip_chord*i/nChord;
                    z_tip[i] = tip_chord*tip.chamber.get(i);
                    x_tip[i] = c_tip*x - s_tip*z_tip[i];
                    z_tip[i] = s_tip*x + c_tip*z_tip[i];
                }

                CartesianMatrix rotate = new CartesianMatrix();
                rotate.getXRotationMatrix(dihedral);

                for(int i = 0; i < nSpan; i++) {
                    int i_outboard = i+1;
                    
                    double le_y_inboard = origin.y + i*dy;
                    double le_y_outboard = origin.y + i_outboard*dy;

                    for(int j = 0; j < nChord; j++) {
                        int j_te = j+1;
                        panels[i][j] = new Panel();
                        double dz_le = (z_tip[j]-z_root[j])/nSpan;
                        double dz_te = (z_tip[j_te]-z_root[j_te])/nSpan;
                        double dx_le = (x_tip[j]-x_root[j])/nSpan;
                        double dx_te = (x_tip[j_te]-x_root[j_te])/nSpan;
                        Cartesian x1y1 = rotate.mult(new Cartesian(x_root[j] + i*dx_le, le_y_inboard, z_root[j]+i*dz_le));
                        Cartesian x2y1 = rotate.mult(new Cartesian(x_root[j_te] + i*dx_te, le_y_inboard, z_root[j_te] + i*dz_te));
                        Cartesian x2y2 = rotate.mult(new Cartesian(x_root[j_te] + i_outboard*dx_te, le_y_outboard, z_root[j_te] + i_outboard*dz_te));
                        Cartesian x1y2 = rotate.mult(new Cartesian(x_root[j] +i_outboard*dx_le, le_y_outboard, z_root[j] + i_outboard*dz_le));
                        panels[i][j].setVertices(x1y1, x2y1, x2y2, x1y2);
                    }
                }
            } else {
                
            }
        }
    }

    /**
     * Add control surface and deflect
     * @param type
     * @param controlYStart
     * @param controlYEnd
     * @param chordFraction
     * @param deflection
     */
    public void addControlSurface(int type, double controlYStart, double controlYEnd, double chordFraction, int deflection) {
        // since constant dy in surface
        double dy = panels[0][0].vertices[2].y-panels[0][0].vertices[0].y;
        double starty = panels[0][0].vertices[0].y+dy/2;
        int[] control = new int[6];
        double test = (controlYStart-starty)/dy;
        control[0] = type;
        control[1] = (int)Math.round((controlYStart-starty)/dy);
        control[2] = (int)Math.round((controlYEnd-starty)/dy);
        if(control[2] >= nSpan) {
            control[2] = nSpan-1;
        }
        control[3] = (int)Math.round((1.0-chordFraction)*nChord);
        if (control[3] < 0) {
            control[3] = 0;
        }
        control[4] = nChord-1;
        // control[5] = deflection;
        controlSurface.add(control);
        deflectControl(controlSurface.size()-1, deflection);
    }
    
    /**
     * Set deflection of contro
     * @param i
     * @param deflection
     */
    public void deflectControl(int i, int deflection) {
        int[] ids = controlSurface.get(i);
        ids[5] = deflection;
        Cartesian LE1 = panels[ids[1]][ids[3]].vertices[0];
        Cartesian LE2 = panels[ids[2]][ids[3]].vertices[0];
        Cartesian axis = LE2.sub(LE1);
        axis.normalize();
        CartesianMatrix M = CartesianMatrix.axisAngle(axis, deflection*MyMath.DEG2RAD);

        for(int j = ids[1]; j <= ids[2]; j++) {
            for(int k = ids[3]; k <= ids[4]; k++) {
                Cartesian[] v = new Cartesian[4];
                Panel p = panels[j][k];
                for(int l = 0; l < 4; l++) {
                    Cartesian r = p.vertices[l].sub(LE1);
                    Cartesian r1 = M.mult(r);
                    v[l] = LE1.add(r1);
                }
                p.setVertices(v[0], v[3], v[2], v[1]);
            }
        }
    }

    public Panel[][] getPanels() {
        return panels;
    }

    public int getChordNumber() {
        return nChord;
    }

    public int getSpanNumber() {
        return nSpan;
    }

    public void setOrigin(Cartesian o) {
        this.origin = o;
    }

    public Cartesian getOrigin() {
        return this.origin;
    }

    public Cartesian getLowerBounds() {
        return this.lowerBounds;
    }

    public Cartesian getUpperBounds() {
        return this.upperBounds;
    }

    public void setName(String s){
        this.name = s;
    }

    public String getName() {
        return name;
    }

    public ArrayList<int[]> getControls() {
        return this.controlSurface;
    }
}