package geometry;

public class GeometryException extends Exception {
    
    public static final int INVALID_INPUT = 1;

    public static final int INVALID_DEFINITION = 2;

    public GeometryException(String message) {
        super(message);
    }
}