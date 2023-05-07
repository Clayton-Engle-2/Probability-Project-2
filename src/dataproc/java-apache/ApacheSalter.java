

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class ApacheSalter {
    private final NormalDistribution normalDistribution = new NormalDistribution();
    private final double CONFIDENCE_LEVEL = 0.99;

    /**
     * Salts an array of x-coordinates by adding a random value to each coordinate using a normal distribution.
     *
     * @param xCoords the array of x-coordinates to salt
     * @return a new array of salted x-coordinates
     */
    public double[] saltCoordinates(double[] xCoords) {
        DescriptiveStatistics stats = new DescriptiveStatistics(xCoords);
        double mean = stats.getMean();
        double std = stats.getStandardDeviation();
        double[] saltedXCoords = new double[xCoords.length];
        for (int i = 0; i < xCoords.length; i++) {
            // Generate a random salt value from a normal distribution with mean and standard deviation of the x-coordinates
            double saltedX = normalDistribution.inverseCumulativeProbability(1 - (1 - CONFIDENCE_LEVEL) / 2) * std + mean;
            saltedXCoords[i] = saltedX;
        }
        return saltedXCoords;
    }
}

