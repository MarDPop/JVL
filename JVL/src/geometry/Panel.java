package geometry;

import utils.Cartesian;

public class Panel {

    public static final double FOUR_PI_INV = 0.25/Math.PI;

    /**
     * Collocation point (1/4 chord along centerline) point where flow needs to be negated aka. Control Point
     */
    private Cartesian collocation; 

    /**
     * Center of panel
     */
    private Cartesian center;

    /**
     * Point located on center of transverse vortex
     */
    private Cartesian transverse; 

    /**
     * Vertices 0 = inboard upstream 1 = inboard downstream 2 = outboard downstream 3 = outboard upstream
     */
    public final Cartesian[] vertices = {new Cartesian(),new Cartesian(),new Cartesian(),new Cartesian()};

    /**
     * 3/4 chord inboard
     */
    private Cartesian A;

    /**
     * 3/4 chord outboard
     */
    private Cartesian B;

    /**
     * transverse vortex segment (A -> B)
     */
    private Cartesian r0;

    /**
     * Normal Vector at collocation point (where needs to be negated)
     */
    private Cartesian normal; 

    /**
     * velocity tangent to panel
     */
    private Cartesian tangentMach;

    /**
     * Normal for compressible flow
     */
    private Cartesian compressibleNormal;

    /**
     * local beta 
     */
    private Cartesian localCompressibleScaling = new Cartesian(1,1,1);;

    /**
     * Circulation Strength
     */
    private double lambda;

    private double area;

    /**
     * Force at ?
     */
    private Cartesian force;

    /**
     * Inducted Velocity at ?
     */
    private Cartesian inducedVelocity;

    public Panel() { }

    /**
     * Simple constructor if not wanting to set vertices
     * @param transverse
     * @param collocation
     * @param n
     * @param area
     */
    public Panel(Cartesian A, Cartesian B, Cartesian n, Cartesian collocation, double lambda) {
        this.transverse = A.add(B).mult(0.5);
        this.A = A;
        this.B = B;
        this.r0 = B.sub(A);
        this.collocation = collocation;
        this.normal = n;
        this.lambda = lambda;
    }

    public Panel(Cartesian x1y1, Cartesian x2y1,Cartesian x2y2,Cartesian x1y2){
        setVertices(x1y1, x2y1, x2y2, x1y2);
    }

    /**
     * Best way to define panel
     * @param x1y1 upstream "inboard" 
     * @param x2y1 downstream "inboard"
     * @param x2y2 downstream "outboard"
     * @param x1y2 upstream "outboard"
     */
    public void setVertices(Cartesian x1y1, Cartesian x2y1,Cartesian x2y2,Cartesian x1y2) {
        this.vertices[0] = x1y1;
        this.vertices[1] = x1y2;
        this.vertices[2] = x2y2;
        this.vertices[3] = x2y1;

        Cartesian edge1 = x2y1.sub(x1y1); // x axis "inboard" edge
        Cartesian edge2 = x1y2.sub(x1y1); // y axis edge
        
        this.normal = edge1.cross(edge2); // z up for positive x and positive y
        this.area = this.normal.getMagnitude(); 
        
        this.normal = this.normal.mult(1/this.area);

        Cartesian midpointLo = x1y1.add(x1y2).mult(0.5);
        Cartesian midpointHi = x2y1.add(x2y2).mult(0.5);

        Cartesian centerline = midpointHi.sub(midpointLo);

        this.collocation = midpointLo.add(centerline.mult(0.75));
        this.center = midpointLo.add(centerline.mult(0.5));
        this.transverse = midpointLo.add(centerline.mult(0.25));

        Cartesian edge3 = x2y2.sub(x1y2); // x axis "outboard" edge
        this.A = this.vertices[0].add(edge1.mult(0.25));
        this.B = this.vertices[1].add(edge3.mult(0.25));
        this.r0 = B.sub(A);
    }

    /**
     * Gets the induced velocity vector using Biot-Savart law at a given point 
     * @param point
     * @param beta
     * @return
     */
    public Cartesian getInducedVelocityFactorAtPointIncompressible(Cartesian point) {

        // Distance to point from Point A (r1) and Point B (r2)
        Cartesian r1 = point.sub(this.A);
        Cartesian r2 = point.sub(this.B);

        // get norm
        Cartesian r1_normalized = r1.mult(1/r1.getMagnitude());
        Cartesian r2_normalized = r2.mult(1/r2.getMagnitude());

        // Contribution from Vortex A-B (transverse segment)
        Cartesian VAB = r1.cross(r2); // temporary storage
        double m = VAB.x*VAB.x+VAB.y*VAB.y+VAB.z*VAB.z;
        if (m < 1e-15){
            VAB.multBy(0.0);
        } else {
            VAB.multBy(r0.dot(r1_normalized.sub(r2_normalized))/m);
        }
        // Contribution from Infinite Vortex Segment A
        Cartesian VA = new Cartesian(0,r1.z,-r1.y);
        m = VA.y*VA.y+VA.z*VA.z;
        if (m < 1e-15){
            VA.multBy(0);
        } else {
            VA.multBy((1+r1_normalized.x)/m);
        }
        // Contribution from Infinite Vortex Segment B
        Cartesian VB = new Cartesian(0,r2.z,-r2.y); 
        m = VB.y*VB.y+VB.z*VB.z;
        if (m < 1e-15){
            VB.multBy(0);
        } else {
            VB.multBy(-(1+r2_normalized.x)/m); // this is inverted from paper (check ) I think this is correct because the B segment vortex vector points away thus dot product is inverted (cross product remains same)
        } 

        return VAB.addTo(VA).addTo(VB).multBy(FOUR_PI_INV);
    }

    /**
     * adds compressiblity effects
     * @param freestream
     * @param speedofsound
     */
    public void addCompressiblity(Cartesian freestream, double speedofsound) {
        this.compressibleNormal = new Cartesian(normal);
        Cartesian edge1 = vertices[3].sub(this.vertices[0]); // x axis "inboard" edge
        Cartesian edge2 = vertices[1].sub(this.vertices[0]);
        for (int i = 0; i < 5; i++){
            tangentMach = freestream.sub(this.compressibleNormal.mult(freestream.dot(this.compressibleNormal)));
            tangentMach.multBy(1/speedofsound);
            localCompressibleScaling.x = 1/(1-tangentMach.x*tangentMach.x);
            localCompressibleScaling.y = 1/(1-tangentMach.y*tangentMach.y);
            localCompressibleScaling.z = 1/(1-tangentMach.z*tangentMach.z);

            this.compressibleNormal = edge1.scale(localCompressibleScaling).cross(edge2.scale(localCompressibleScaling));
            this.compressibleNormal.normalize();
        }

    }

        /**
     * Gets the induced velocity vector using Biot-Savart law at a given point and mach factor (1-Mach^2)
     * @param point
     * @param beta
     * @return
     */
    public Cartesian getInducedVelocityFactorAtPointCompressible(Cartesian point, Cartesian Scaling) {

        // Distance to point from Point A (r1) and Point B (r2)
        Cartesian r1 = point.sub(this.A);
        Cartesian r2 = point.sub(this.B);

        r1.scaleBy(localCompressibleScaling);
        r2.scaleBy(localCompressibleScaling);

        // get norm
        Cartesian r1_normalized = r1.mult(1/r1.getMagnitude());
        Cartesian r2_normalized = r2.mult(1/r2.getMagnitude());

        // Contribution from Vortex A-B (transverse segment)
        Cartesian VAB = r1.cross(r2); // temporary storage
        double m = VAB.x*VAB.x+VAB.y*VAB.y+VAB.z*VAB.z;
        if (m < 1e-15){
            VAB.multBy(0.0);
        } else {
            VAB.multBy(r0.dot(r1_normalized.sub(r2_normalized))/m);
        }
        // Contribution from Infinite Vortex Segment A
        Cartesian VA = new Cartesian(0,r1.z,-r1.y);
        m = VA.y*VA.y+VA.z*VA.z;
        if (m < 1e-15){
            VA.multBy(0);
        } else {
            VA.multBy((1+r1_normalized.x)/m);
        }
        // Contribution from Infinite Vortex Segment B
        Cartesian VB = new Cartesian(0,r2.z,-r2.y); 
        m = VB.y*VB.y+VB.z*VB.z;
        if (m < 1e-15){
            VB.multBy(0);
        } else {
            VB.multBy(-(1+r2_normalized.x)/m); // this is inverted from paper (check ) I think this is correct because the B segment vortex vector points away thus dot product is inverted (cross product remains same)
        } 

        return VAB.addTo(VA).addTo(VB).multBy(FOUR_PI_INV).scaleBy(localCompressibleScaling);
    }

    /* Getters and Setters */

    public void setForce(Cartesian force) {
        this.force = force;
    }

    public Cartesian getForce() {
        return this.force;
    }

    public void setInducedVelocity(Cartesian inducedVelocity) {
        this.inducedVelocity = inducedVelocity;
    }

    public Cartesian getInducedVelocity() {
        return this.inducedVelocity;
    }

    // these setters shouldn't be used, use the set vertices

    public Cartesian getCollocationPoint() {
        return this.collocation;
    }

    public void setCollocationPoint(double x, double y, double z) {
        this.collocation.set(x,y,z);
    }

    public Cartesian getA() {
        return A;
    }

    public Cartesian getB() {
        return B;
    }

    public Cartesian getVortexVector() {
        return r0;
    }

    public Cartesian getTransversePoint() {
        return this.transverse;
    }

    public void setTransversePoint(double x, double y, double z) {
        this.transverse.set(x,y,z);
    }

    public void setNormal(double x, double y, double z) {
        this.normal.set(x,y,z);
    }

    public void setNormal(Cartesian n) {
        this.normal = n;
    }
    
    public Cartesian getNormal() {
        return this.normal;
    }
    
    public Cartesian getCompressibleNormal() {
        return this.compressibleNormal;
    }

    public void setCirculation(double l) {
        this.lambda = l;
    }

    public double getCirculation() {
        return this.lambda;
    }

    public void setArea(double area) {
        this.area = area;
    }

    public double getArea() {
        return this.area;
    }

    public Cartesian getCenter() {
        return this.center;
    }

}