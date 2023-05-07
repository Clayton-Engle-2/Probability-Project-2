package research.optimize;

import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class ApachePlotter {
    private double xStart;
    private double xEnd;
    private double frequency;

    public ApachePlotter(double xStart, double xEnd, double frequency) {
        this.xStart = xStart;
        this.xEnd = xEnd;
        this.frequency = frequency;
    }

    /**
     * Generates an array of x and y values for the function y = sin(x).
     *
     * @return an array of (x,y) pairs for y = sin(x)
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
     * Creates a new CSV file and writes the given data to it.
     *
     * @param filename the name of the CSV file to create
     * @param data     the data to write to the file
     * @throws IOException if there is an error creating or writing to the file
     */
    public void writeToFile(String filename, double[][] data) throws IOException {
        FileWriter writer = new FileWriter(filename);
        CSVPrinter printer = new CSVPrinter(writer, CSVFormat.DEFAULT.withHeader("x", "y"));

        for (int i = 0; i < data.length; i++) {
            printer.printRecord(data[i][0], data[i][1]);
        }

        printer.close();
        writer.close();
    }

    /**
     * Generates an array of (x,y) pairs for the function y = sin(x) and writes the data to a CSV file.
     *
     * @param filename the name of the CSV file to create
     * @throws IOException if there is an error creating or writing to the file
     */
    public void newSinPlot(String filename) throws IOException {
        double[][] data = sinx();
        writeToFile(filename, data);
    }

    /**
     * Generates an array of (x,y) pairs for the function y = f(x) using cubic spline interpolation.
     *
     * @param function the function to interpolate
     * @param numPoints the number of points to generate
     * @return an array of (x,y) pairs for y = f(x)
     */
    public double[][] interpolate(PolynomialSplineFunction function, int numPoints) {
        double[][] data = new double[numPoints][2];
        double step = (xEnd - xStart) / (numPoints - 1);

        for (int i = 0; i < numPoints; i++) {
            double x = xStart + i * step;
            double y = function.value(x);
            data[i][0] = x;
            data[i][1] = y;
        }

        return data;
    }

    /**
     * Generates an array of (x,y) pairs for the function y = f(x) using cubic spline interpolation and writes the data to a CSV file.
     *
          * @param filename the name of the CSV file to create
     * @param function the function to interpolate
     * @param numPoints the number of points to generate
     * @throws IOException if there is an error creating or writing to the file
     */
    public void newInterpolationPlot(String filename, PolynomialSplineFunction function, int numPoints) throws IOException {
        double[][] data = interpolate(function, numPoints);
        writeToFile(filename, data);
    }

    /**
     * Generates a polynomial spline function that interpolates the given data points.
     *
     * @param data the data points to interpolate
     * @return a polynomial spline function that interpolates the data points
     */
    public PolynomialSplineFunction interpolate(double[][] data) {
        SplineInterpolator interpolator = new SplineInterpolator();
        double[] x = new double[data.length];
        double[] y = new double[data.length];

        for (int i = 0; i < data.length; i++) {
            x[i] = data[i][0];
            y[i] = data[i][1];
        }

        return interpolator.interpolate(x, y);
    }
}


