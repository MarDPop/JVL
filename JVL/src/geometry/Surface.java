package geometry;

import utils.Cartesian;

public class Surface {

    String name;
    
    Panel[][] panels;

    int nChord;

    int nSpan;

    Cartesian origin = new Cartesian();

    Cartesian lowerBounds = new Cartesian();

    Cartesian upperBounds = new Cartesian();

    public Surface() {}

    public Surface(Cartesian origin, double root_chord, double span,int nChord, int nSpan, double angle, double taper, double sweep) {

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
        double spanRatioTaper = Math.tan(taper)*d_span;
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

    public Surface(Cartesian LE1, Curve chord1, Cartesian LE2, Curve chord2, int nSpan) {

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
}