package geometry;

public class GeometryException extends Exception {
    
    /**
     *
     */
    private static final long serialVersionUID = -6094029487515286972L;

    public static final int INVALID_INPUT = 1;

    public static final int INVALID_DEFINITION = 2;

    public GeometryException(String message) {
        super(message);
    }
}