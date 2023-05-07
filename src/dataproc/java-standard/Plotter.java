
import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;

public class Plotter {

    /**
     * Constructs a new Plotter object with the given x start, x end,
     * and frequency.
     *
     * @param xStart    the starting value for x
     * @param xEnd      the ending value for x
     * @param frequency the frequency of the sine wave
     */
    public Plotter(double xStart, double xEnd, double frequency) {
        this.xStart = xStart;
        this.xEnd = xEnd;
        this.frequency = frequency;
    }

    /**
     * Generates an array of data points for a sine wave.
     *
     * @return a 2D array of doubles where each row represents an (x, y) pair
     */
    public double[][] sinx() {
        int numPoints = (int) Math.ceil((xEnd - xStart) / frequency) + 1;
        double[][] data = new double[numPoints][2];

        for (int i = 0; i < numPoints; i++) {
            double x = xStart + i * frequency;
            double y = Math.sin(x);
            data[i][0] = x;
            data[i][1] = y;
        }

        return data;
    }

    /**
     * Writes a 2D array of data points to a CSV file with the given filename.
     *
     * @param filename the name of the file to write to
     * @param data     the data points to write to the file
     * @throws IOException if an I/O error occurs while writing to the file
     */
    public void writeToFile(String filename, double[][] data) throws IOException {
        try (FileWriter writer = new FileWriter(filename)) {
            writer.write("x,y\n");
            for (double[] row : data) {
                writer.write(String.format(Locale.US, "%.8f,%.8f%n", row[0], row[1]));
            }
        }
    }

    /**
     * Generates a sine wave and writes the data points to a CSV file with
     * the given filename.
     *
     * @param filename the name of the file to write to
     * @throws IOException if an I/O error occurs while writing to the file
     */
    public void newSinPlot(String filename) throws IOException {
        double[][] data = sinx();
        writeToFile(filename, data);
    }
}
