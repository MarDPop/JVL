package geometry;

import java.util.ArrayList;
import java.util.HashMap;

public class Airfoil {
    ArrayList<Double> thickness = new ArrayList<>();

    ArrayList<Double> chamber = new ArrayList<>();

    ArrayList<Double> upperSurfaceY = new ArrayList<>();

    ArrayList<Double> upperSurfaceX = new ArrayList<>();
    
    ArrayList<Double> lowerSurfaceX = new ArrayList<>();

    ArrayList<Double> lowerSurfaceY = new ArrayList<>();

    private static final HashMap<String,double[]> table = new HashMap<>();

    static {
        // airfoil thickness, position of max chord, k1, k2/k1 (-1 if not reflex)
        table.put("210",new double[]{0.05,0.0580,361.400,-1});
        table.put("220",new double[]{0.10,0.1260,51.640,-1});
        table.put("230",new double[]{0.15,0.2025,15.957,-1});
        table.put("240",new double[]{0.20,0.2900,6.643,-1});
        table.put("250",new double[]{0.25,0.3910,3.230,-1});

        table.put("221",new double[]{0.10,0.1300,51.990,0.000764});
        table.put("231",new double[]{0.15,0.2170,15.793,0.00677});
        table.put("241",new double[]{0.20,0.3180,6.520,0.0303});
        table.put("251",new double[]{0.25,0.4410,3.191,0.1355});
    }

    public Airfoil() {
    }

    public Airfoil(String Naca, int divisions) {
        if(Naca.length() == 4) {
            int m = Naca.charAt(0)-48;
            int p = Naca.charAt(1)-48;
            int thick = Integer.parseInt(Naca.substring(2));
            setNACA4Series(thick/100.0, m/100.0, p/10.0, divisions);
        } else {
            if(Naca.length() == 5) {
                double[] vals = table.get(Naca.substring(0,3));
                int thick = Integer.parseInt(Naca.substring(3));
                setNACA5Series(thick/100.0, vals[0], vals[1], vals[2], vals[3], divisions);
            } else {
                
            }
        }
    }
    
    public void setNACA4Series(double maxThickness, double maxChamber, double chamberPosition, int divisions) {
        double dx = 1.0/divisions;
        chamber.add(0.0);
        thickness.add(0.0);
        for(int i = 1; i < divisions; i++) {
            double x = i*dx;
            thickness.add(maxThickness/0.2*(0.2969*Math.sqrt(x)+x*(-0.126+x*(-0.3516+x*(0.2843-x*0.1015)))));
            double theta;
            if(x < chamberPosition) {
                double xp = x/chamberPosition;
                chamber.add(maxChamber*(xp*(2-xp)));
                theta = Math.atan(maxChamber*2*(1-xp)/(chamberPosition));
                
            } else {
                double p2 = 1-chamberPosition;
                chamber.add(maxChamber*(1 - 2*chamberPosition*(1 - x) - x*x)/(p2*p2));
                theta = Math.atan(maxChamber*2*(chamberPosition-x)/(p2*p2));
            }
            upperSurfaceX.add(x-thickness.get(i)*Math.sin(theta));
            upperSurfaceY.add(chamber.get(i)+thickness.get(i)*Math.cos(theta));
            upperSurfaceX.add(x+thickness.get(i)*Math.sin(theta));
            upperSurfaceY.add(chamber.get(i)-thickness.get(i)*Math.cos(theta));
        }
        chamber.add(0.0);
        thickness.add(0.0);
    }

    public void setNACA5Series(double maxThickness, double p, double m, double k1, double k2k1, int divisions) {
        double dx = 1.0/divisions;
        chamber.add(0.0);
        thickness.add(0.0);
        if (k2k1 < 0) {
            for(int i = 1; i < divisions; i++) {
                double x = i*dx;
                thickness.add(maxThickness/0.2*(0.2969*Math.sqrt(x)+x*(-0.126+x*(-0.3516+x*(0.2843-x*0.1015)))));
                double theta;
                if(x < p) {
                    chamber.add( k1/6*x*( m*m*(3-m) + x*(-3*m + x) ) );
                    theta = Math.atan(k1/6*(m*m*(3-m)+2*x*(3*m-x/3)));

                } else {
                    chamber.add(k1/6*m*m*m*(1-x));
                    theta = Math.atan(-k1/6*m*m*m);
                }
                upperSurfaceX.add(x-thickness.get(i)*Math.sin(theta));
                upperSurfaceY.add(chamber.get(i)+thickness.get(i)*Math.cos(theta));
                upperSurfaceX.add(x+thickness.get(i)*Math.sin(theta));
                upperSurfaceY.add(chamber.get(i)-thickness.get(i)*Math.cos(theta));
            }
        } else {
            for(int i = 1; i < divisions; i++) {
                double x = i*dx;
                thickness.add(maxThickness/0.2*(0.2969*Math.sqrt(x)+x*(-0.126+x*(-0.3516+x*(0.2843-x*0.1015)))));
                double theta;
                if(x < p) {
                    double r = x-p;
                    double r2 = 1-p;
                    chamber.add( k1/6*(r*r*r + p*p*p*(1-x) - k2k1*r2*r2*r2*x)  );
                    theta = Math.atan( k1/6*(k2k1*r2*r2*r2 - p*p*p + r*r*3) );

                } else {
                    double r = x-p;
                    double r2 = 1-p;
                    chamber.add( k1/6*(p*p*p*(1-x) - k2k1*(r2*r2*r2*x-r*r*r))  );
                    theta = Math.atan( k1/6*(k2k1*(-r2*r2*r2 + r*r*3) - p*p*p ) );
                }
                upperSurfaceX.add(x-thickness.get(i)*Math.sin(theta));
                upperSurfaceY.add(chamber.get(i)+thickness.get(i)*Math.cos(theta));
                upperSurfaceX.add(x+thickness.get(i)*Math.sin(theta));
                upperSurfaceY.add(chamber.get(i)-thickness.get(i)*Math.cos(theta));
            }
        }
        chamber.add(0.0);
        thickness.add(0.0);
    }


    public ArrayList<Double> getChamber() {
        return chamber;
    }

}