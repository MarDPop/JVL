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

    public Surface(Cartesian origin, double chord, double span,int nChord, int nSpan) {

        if(nChord < 0 || nSpan < 0 ) {
            nChord = Math.abs(nChord);
            nSpan = Math.abs(nSpan);
            System.out.println("Input number of divisions must be positive, absolute value used.");
        }

        this.nChord = nChord;
        this.nSpan = nSpan;

        double dx = chord/nChord;
        double dy = span/nSpan;

        panels = new Panel[nSpan][nChord];
        
        for(int i = 0; i < nSpan; i++) {
            for(int j = 0; j < nChord; j++) {
                panels[i][j] = new Panel();
                Cartesian x1y1 = new Cartesian(origin.x + j*dx,origin.y + i*dy,0);
                Cartesian x2y1 = new Cartesian(origin.x +(j+1)*dx,origin.y +i*dy,0);
                Cartesian x2y2 = new Cartesian(origin.x +(j+1)*dx,origin.y +(i+1)*dy,0);
                Cartesian x1y2 = new Cartesian(origin.x +j*dx,origin.y +(i+1)*dy,0);
                panels[i][j].setVertices(x1y1, x2y1, x2y2, x1y2);
            }
        }

        this.origin = origin;

        if(chord > 0) {
            this.lowerBounds.x = origin.x;
            this.upperBounds.x = origin.x+chord;
        } else {
            this.lowerBounds.x = origin.x+chord;
            this.upperBounds.x = origin.x;
        }
        
        if(span > 0) {
            this.lowerBounds.y = origin.y;
            this.upperBounds.y = origin.y+span;
        } else {
            this.lowerBounds.y = origin.y+span;
            this.upperBounds.y = origin.y;
        }

        this.lowerBounds.z = origin.z;
        this.upperBounds.z = origin.z;

    }

    public Surface(Cartesian origin, double span, double chord, int nSpan, int nChord, double taper, double sweep, double dihedral) {

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