package app;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JPanel;
import javax.swing.SwingConstants;

import javax.swing.JLabel;
import java.awt.Graphics;
import java.awt.Color;
import java.awt.Font;

public class PlotPanel extends JPanel {

    ArrayList<ArrayList<Double>> x = new ArrayList<>();

    ArrayList<ArrayList<Double>> y = new ArrayList<>();

    ArrayList<ArrayList<Double>> z = new ArrayList<>();

    ArrayList<String> legends;

    double max_x,min_x,max_y,min_y,max_z,min_z;

    Map<String,Object> options;

    public static final Map<String,Object> DEFAULT_OPTIONS = new HashMap<>();

    private static final int DEFAULT_MARGINS = 100;

    boolean is2D;

    Graphics g;

    static {
        DEFAULT_OPTIONS.put("Color",Color.black);
        DEFAULT_OPTIONS.put("BackgroundColor",Color.white);
        DEFAULT_OPTIONS.put("MajorTickNumber",8);
    }

    public PlotPanel() {}

    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        this.g = g;
        g.setColor(Color.white);
        g.fillRect(DEFAULT_MARGINS,DEFAULT_MARGINS,this.getWidth()-2*DEFAULT_MARGINS,this.getHeight()-2*DEFAULT_MARGINS);
        replot();
    }

    /**
     * Simple plotter
     * @param x
     * @param y
     */
    public void plot(ArrayList<Double> x, ArrayList<Double> y, Map<String,Object> options) throws IllegalArgumentException {

        if(x.size() != y.size()) {
            throw new IllegalArgumentException("size of arrays must match");
        }

        this.x.add(x);
        this.y.add(y);
        min_x = x.get(0);
        max_x = x.get(0);
        min_y = y.get(0);
        max_y = y.get(0);
        this.options = options;
        for(double d : x) {
            if(d < min_x) {
                min_x = d;
            } 
            if(d > max_x) {
                max_x = d;
            }
        }
        for(double d : y) {
            if(d < min_y) {
                min_y = d;
            } 
            if(d > max_y) {
                max_y = d;
            }
        }
        is2D = true;

        repaint();
    }

    private void replot() {

        this.removeAll();

        int pixelsHeight = getHeight() - 2*DEFAULT_MARGINS;
        int pixelsWidth = getWidth() - 2*DEFAULT_MARGINS;
        int bottomY = getHeight()-DEFAULT_MARGINS;
        int rightX = getWidth() - DEFAULT_MARGINS;
        
        g.setColor(Color.black);

        double xRange = max_x - min_x;
        double yRange = max_y - min_y;

        g.drawLine(DEFAULT_MARGINS,DEFAULT_MARGINS,DEFAULT_MARGINS,bottomY);
        g.drawLine(DEFAULT_MARGINS,bottomY,rightX,bottomY);
        
        int tmp1 = bottomY + 10;
        int tmp2 = DEFAULT_MARGINS-10;
        int divisions = (int)DEFAULT_OPTIONS.get("MajorTickNumber");
        // g.setFont(new Font(Font.SANS_SERIF,Font.BOLD,10));
        for(int i = 0; i <= divisions; i++) {
            // y axis
            int p = DEFAULT_MARGINS+(pixelsHeight*i)/divisions;
            g.drawLine(tmp2,p,DEFAULT_MARGINS,p);
            double y = max_y - i*yRange/divisions;
            JLabel l = new JLabel(Double.toString(y),SwingConstants.RIGHT);
            l.setSize(50,10);
            l.setLocation(tmp2-60,p-5);
            this.add(l);
            // x axis
            p = DEFAULT_MARGINS+(pixelsWidth*i)/divisions;
            g.drawLine(p,tmp1,p,bottomY);
            double x = min_x + i*xRange/divisions;
            l = new JLabel(Double.toString(x),SwingConstants.CENTER);
            l.setSize(50,10);
            l.setLocation(p-25,tmp1+10);
            this.add(l);
        }

        if(options.containsKey("Title")) {
            JLabel t = new JLabel((String)options.get("Title"),SwingConstants.CENTER);
            t.setLocation(getWidth()/2-200,30);
            t.setSize(400,40);
            t.setFont(new Font(Font.SANS_SERIF,Font.BOLD,40));
            this.add(t);
        }

        if(options.containsKey("Color")) {
            g.setColor((Color)options.get("Color"));
        } else {
            g.setColor(Color.blue);
        }

        for(int i = 0; i < x.size(); i++) {
            ArrayList<Double> xs = x.get(i);
            ArrayList<Double> ys = y.get(i);
            
            if(xs.size() > 1) {
                for(int j = 1; j < xs.size(); j++) {
                    g.drawLine((int)(DEFAULT_MARGINS+pixelsWidth*((xs.get(j-1)-min_x)/xRange)),(int)(bottomY-pixelsHeight*((ys.get(j-1)-min_y)/yRange)),(int)(DEFAULT_MARGINS+pixelsWidth*((xs.get(j)-min_x)/xRange)),(int)(bottomY-pixelsHeight*((ys.get(j)-min_y)/yRange)));
                }
            }
            
        }
      
    }

    
}