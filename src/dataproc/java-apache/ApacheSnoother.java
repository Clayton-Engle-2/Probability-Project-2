package apache;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class ApacheSmoother {

	public static double[] smooth(double[] x, double[] y, int windowSize, int polynomialOrder) {
		// Check that both x and y have the same length
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y must have the same length");
		}

		// Validate windowSize and polynomialOrder
		if (windowSize % 2 == 0) {
			throw new IllegalArgumentException("windowSize must be odd");
		}
		if (windowSize < polynomialOrder + 1) {
			throw new IllegalArgumentException("windowSize must be at least polynomialOrder + 1");
		}

		// Create an array to store the smoothed y values
		double[] smoothedY = new double[y.length];

		// Calculate the Savitzky-Golay coefficients
		double[] coefficients = savitzkyGolayCoefficients(windowSize, polynomialOrder);

		// Apply the Savitzky-Golay filter
		int halfWindowSize = windowSize / 2;
		for (int i = halfWindowSize; i < y.length - halfWindowSize; i++) {
			double sum = 0;
			for (int j = 0; j < windowSize; j++) {
				sum += y[i - halfWindowSize + j] * coefficients[j];
			}
			smoothedY[i] = sum;
		}

		// Fill in the edges of the data with the original y values
		System.arraycopy(y, 0, smoothedY, 0, halfWindowSize);
		System.arraycopy(y, y.length - halfWindowSize, smoothedY, y.length - halfWindowSize, halfWindowSize);

		// Return the smoothed y values
		return smoothedY;
	}

	public static double[] savitzkyGolayCoefficients(int windowSize, int polynomialOrder) {
		// Create the Vandermonde matrix
		RealMatrix vandermondeMatrix = MatrixUtils.createRealMatrix(windowSize, polynomialOrder + 1);
		for (int i = 0; i < windowSize; i++) {
			for (int j = 0; j <= polynomialOrder; j++) {
				vandermondeMatrix.setEntry(i, j, Math.pow(i - windowSize / 2, j));
			}
		}

		// Compute the pseudoinverse of the Vandermonde matrix
		RealMatrix pseudoinverse = MatrixUtils.inverse(vandermondeMatrix.transpose().multiply(vandermondeMatrix))
				.multiply(vandermondeMatrix.transpose());

		// Return the first row of the pseudoinverse as the Savitzky-Golay coefficients
		double[] coefficients = new double[windowSize];
		for (int i = 0; i < windowSize; i++) {
			coefficients[i] = pseudoinverse.getEntry(0, i);
		}

		return coefficients;
	}

}
