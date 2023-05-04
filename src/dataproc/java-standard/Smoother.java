package standard;

public class Smoother {
	
	public static double[] smooth(double[] x, double[] y, double alpha) {
	    // Check that both x and y have the same length
	    if (x.length != y.length) {
	        throw new IllegalArgumentException("x and y must have the same length");
	    }

	    // Create an array to store the smoothed y values
	    double[] smoothedY = new double[y.length];

	    // Initialize the smoothedY with the first y value
	    smoothedY[0] = y[0];

	    // Apply the exponentially weighted moving average to the rest of the data
	    for (int i = 1; i < y.length; i++) {
	        smoothedY[i] = alpha * y[i] + (1 - alpha) * smoothedY[i - 1];
	    }

	    // Return the smoothed y values
	    return smoothedY;
	}


}
