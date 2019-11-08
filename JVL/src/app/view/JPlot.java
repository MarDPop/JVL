package app.view;

import java.util.ArrayList;
import java.util.Map;

import javax.swing.JFrame;

public class JPlot extends JFrame {

    /**
     *
     */
    private static final long serialVersionUID = -8942536461877043607L;
    
    PlotPanel plotArea;

    public JPlot() {
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.setSize(800,800);

        this.plotArea = new PlotPanel();

        getContentPane().add(plotArea);

        setVisible(true);
    }

    public JPlot(ArrayList<Double> x, ArrayList<Double> y, Map<String,Object> options) {
        this();
        plotArea.plot(x,y,options);
    }

    

}