package app;

public class Options {
    public boolean compressible;

    public boolean scaleCompressible;

    public boolean unsteady;

    public Options() {
        setDefaults();
    }

    public void setDefaults() {
        this.compressible = false;
        this.scaleCompressible = true;
        this.unsteady = false;
    }

}