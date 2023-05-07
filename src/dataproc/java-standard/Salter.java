

import java.util.Random;

public class Salter {
    private final Random random = new Random();

    /**
     * Salts an array of x-coordinates by adding a random value to each coordinate within the specified range.
     *
     * @param xCoords the array of x-coordinates to salt
     * @param range   the maximum distance from the original value that the salted value can be
     * @return a new array of salted x-coordinates
     */
    public double[] saltCoordinates(double[] xCoords, double range) {
        double[] saltedXCoords = new double[xCoords.length];
        for (int i = 0; i < xCoords.length; i++) {
            // Calculate a random salt value within the specified range and add it to the original value
            double saltedX = xCoords[i] + (random.nextDouble() * range * 2) - range;
            saltedXCoords[i] = saltedX;
        }
        return saltedXCoords;
    }
}

