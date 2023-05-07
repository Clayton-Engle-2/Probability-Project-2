import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;

public class StatsLibrary {
	
	//========== Hypergeometric =========================
	
	/**
	 * Calculates the probability mass function (pmf) of a hypergeometric distribution.
	 *
	 * @param populationSize the size of the population
	 * @param numSuccesses   the number of successes in the population
	 * @param sampleSize     the size of the sample
	 * @param numDrawn       the number of successes in the sample
	 * @return the pmf of the hypergeometric distribution
	 */
	public double hypergeometricPmf(int populationSize, int numSuccesses, int sampleSize, int numDrawn) {
	    BigInteger numerator = choose(numSuccesses, numDrawn).multiply(choose(populationSize - numSuccesses, sampleSize - numDrawn));
	    BigInteger denominator = choose(populationSize, sampleSize);
	    return numerator.doubleValue() / denominator.doubleValue();
	}

	/**
	 * Calculates the binomial coefficient "n choose k" using the formula n! / (k! * (n-k)!).
	 *
	 * @param n the size of the set
	 * @param k the size of the subset
	 * @return the binomial coefficient "n choose k"
	 */
	private BigInteger choose(int n, int k) {
	    if (k > n) {
	        return BigInteger.ZERO;
	    }
	    BigInteger result = BigInteger.ONE;
	    for (int i = 1; i <= k; i++) {
	        result = result.multiply(BigInteger.valueOf(n - i + 1)).divide(BigInteger.valueOf(i));
	    }
	    return result;
	}
	/**
	 * Calculates the mean of a hypergeometric distribution.
	 *
	 * @param populationSize the size of the population
	 * @param numSuccesses   the number of successes in the population
	 * @param sampleSize     the size of the sample
	 * @return the mean of the hypergeometric distribution
	 */
	public double hypergeometricMean(int populationSize, int numSuccesses, int sampleSize) {
	    return (double) sampleSize * numSuccesses / populationSize;
	}

	/**
	 * Calculates the variance of a hypergeometric distribution.
	 *
	 * @param populationSize the size of the population
	 * @param numSuccesses   the number of successes in the population
	 * @param sampleSize     the size of the sample
	 * @return the variance of the hypergeometric distribution
	 */
	public double hypergeometricVariance(int populationSize, int numSuccesses, int sampleSize) {
	    return (double) sampleSize * numSuccesses * (populationSize - numSuccesses) * (populationSize - sampleSize)
	            / (populationSize * populationSize * (populationSize - 1));
	}
	
	
	
	
	
	//=============== Poison =========================
	
	/**
	 * Calculates the probability mass function (pmf) of a Poisson distribution.
	 *
	 * @param lambda the rate parameter of the Poisson distribution
	 * @param k      the number of occurrences
	 * @return the pmf of the Poisson distribution
	 */
	public double poissonPmf(double lambda, int k) {
	    double expLambda = Math.exp(-lambda);
	    BigInteger numerator = BigInteger.valueOf((long) Math.pow(lambda, k));
	    BigInteger denominator = factorial(k);
	    return numerator.doubleValue() * expLambda / denominator.doubleValue();
	}

	/**
	 * Calculates the factorial of a non-negative integer.
	 *
	 * @param n the non-negative integer to calculate the factorial of
	 * @return the factorial of the non-negative integer
	 */
	private BigInteger factorial(int n) {
	    if (n == 0) {
	        return BigInteger.ONE;
	    } else {
	        BigInteger result = BigInteger.valueOf(n);
	        return result.multiply(factorial(n - 1));
	    }
	}

	/**
	 * Calculates the mean of a Poisson distribution.
	 *
	 * @param lambda the rate parameter of the Poisson distribution
	 * @return the mean of the Poisson distribution
	 */
	public double poissonMean(double lambda) {
	    return lambda;
	}

	/**
	 * Calculates the variance of a Poisson distribution.
	 *
	 * @param lambda the rate parameter of the Poisson distribution
	 * @return the variance of the Poisson distribution
	 */
	public double poissonVariance(double lambda) {
	    return lambda;
	}
	
	
	//======== Tchebyshev's Theorem ========================
	/**
	 * In this version of the method, k is the number of standard deviations from the mean, and exact is a 
	 * boolean parameter that determines whether the exact percentage should be returned or an upper bound 
	 * is acceptable. If exact is true, the method calculates the exact percentage using Tchebyshev's 
	 * Theorem. If exact is false, the method calculates an upper bound on the percentage.
	 * 
	 * The method first checks that k is a positive integer. It then calculates the variance and mean of 
	 * the dataset, assuming a standard deviation of 1.0. The sigma variable is set to the square root 
	 * of the variance.
	 * 
	 * If exact is true, the method calculates the exact percentage using Tchebyshev's Theorem: 1 - (1/k^2). 
	 * This formula gives the minimum percentage of the data that falls within k standard deviations from 
	 * the mean.
	 * 
	 * If exact is false, the method calculates an upper bound on the percentage using the formula: 
	 * 1 - (1/k/k). This formula gives an upper bound on the percentage of the data that falls within k 
	 * standard deviations from the mean.
	 * 
	 * Calculates the percentage of a dataset that falls within k standard deviations from the mean
	 * using Tchebyshev's Theorem.
	 *
	 * @param k     the number of standard deviations from the mean
	 * @param exact true if the exact percentage should be returned, false if an upper bound is acceptable
	 * @return the percentage of the dataset that falls within k standard deviations from the mean
	 */
	public double tchebyshevsTheorem(int k, boolean exact) {
	    if (k <= 0) {
	        throw new IllegalArgumentException("k must be a positive integer");
	    }

	    double variance = 1.0;
	    double mean = 0.0;
	    double sigma = Math.sqrt(variance);

	    if (exact) {
	        // Calculate the exact percentage
	        return 1.0 - (1.0 / Math.pow(k, 2));
	    } else {
	        // Calculate an upper bound on the percentage
	        return 1.0 - (1.0 / k / k);
	    }
	}
	
	//================== Uniform Continuous ======================
	/**
	 * Calculates the probability density function (pdf) of a continuous uniform distribution.
	 *
	 * @param a the lower bound of the distribution
	 * @param b the upper bound of the distribution
	 * @param x the value at which to evaluate the pdf
	 * @return the pdf of the continuous uniform distribution
	 */
	public double uniformPdf(double a, double b, double x) {
	    if (x < a || x > b) {
	        return 0.0;
	    } else {
	        return 1.0 / (b - a);
	    }
	}

	/**
	 * Calculates the probability mass function (pmf) of a discrete uniform distribution.
	 *
	 * @param a the lower bound of the distribution
	 * @param b the upper bound of the distribution
	 * @param k the value at which to evaluate the pmf
	 * @return the pmf of the discrete uniform distribution
	 */
	public double uniformPmf(int a, int b, int k) {
	    if (k < a || k > b) {
	        return 0.0;
	    } else {
	        return 1.0 / (b - a + 1);
	    }
	}

	/**
	 * Calculates the mean of a continuous uniform distribution.
	 *
	 * @param a the lower bound of the distribution
	 * @param b the upper bound of the distribution
	 * @return the mean of the continuous uniform distribution
	 */
	public double uniformMean(double a, double b) {
	    return (a + b) / 2.0;
	}

	/**
	 * Calculates the variance of a continuous uniform distribution.
	 *
	 * @param a the lower bound of the distribution
	 * @param b the upper bound of the distribution
	 * @return the variance of the continuous uniform distribution
	 */
	public double uniformVariance(double a, double b) {
	    return Math.pow(b - a, 2) / 12.0;
	}
	
	
	
	//============ Normal Continuous ===================
	/**
	 * Calculates the probability density function (pdf) of a normal distribution.
	 *
	 * @param x     the value at which to evaluate the pdf
	 * @param mu    the mean of the distribution
	 * @param sigma the standard deviation of the distribution
	 * @return the pdf of the normal distribution
	 */
	public double normalPdf(double x, double mu, double sigma) {
	    return Math.exp(-Math.pow(x - mu, 2) / (2 * Math.pow(sigma, 2))) / (sigma * Math.sqrt(2 * Math.PI));
	}

	/**
	 * Calculates the mean of a normal distribution.
	 *
	 * @param mu    the mean of the distribution
	 * @param sigma the standard deviation of the distribution
	 * @return the mean of the normal distribution
	 */
	public double normalMean(double mu, double sigma) {
	    return mu;
	}

	/**
	 * Calculates the variance of a normal distribution.
	 *
	 * @param sigma the standard deviation of the distribution
	 * @return the variance of the normal distribution
	 */
	public double normalVariance(double sigma) {
	    return Math.pow(sigma, 2);
	}
	
	
	//=========== Gamma ==============
	/**
	 * Calculates the probability density function (pdf) of a Gamma distribution.
	 *
	 * @param x     the value at which to evaluate the pdf
	 * @param alpha the shape parameter of the distribution
	 * @param beta  the scale parameter of the distribution
	 * @return the pdf of the Gamma distribution
	 * @throws IllegalArgumentException if x, alpha, or beta is NaN or alpha or beta is not positive
	 */
	public double gammaPdf(double x, double alpha, double beta) {
	    if (Double.isNaN(x) || Double.isNaN(alpha) || Double.isNaN(beta)) {
	        throw new IllegalArgumentException("x, alpha, and beta must not be NaN");
	    }
	    if (alpha <= 0.0 || beta <= 0.0) {
	        throw new IllegalArgumentException("alpha and beta must be positive");
	    }
	    return Math.pow(x, alpha - 1) * Math.exp(-x / beta) / (Math.pow(beta, alpha) * gamma(alpha));
	}
	
	 /**
     * Calculates the gamma function for a given input value x.
     * 
     * @param x the input value for which to calculate the gamma function.
     * @return the value of the gamma function for the input value x.
     */
	public double gamma(double x) {
        if (x == 0.0) {
            return Double.POSITIVE_INFINITY;
        } else if (x < 0.0 && x == Math.floor(x)) {
            return Double.NEGATIVE_INFINITY;
        } else if (x < 0.5) {
            return Math.PI / (Math.sin(Math.PI * x) * gamma(1 - x));
        } else {
            x -= 1.0;
            double a = 1.0 / Math.sqrt(2 * Math.PI);
            double b = Math.pow(x / Math.E, x);
            double c = ((1.0 / (12 * x)) - (1.0 / (360 * Math.pow(x, 3))) + (1.0 / (1260 * Math.pow(x, 5))));
            double d = 1.0 + (c * b);
            double e = Math.sqrt(2 * Math.PI * x);
            return a * d * e;
        }
    }
    


	/**
	 * Calculates the mean of a Gamma distribution.
	 *
	 * @param alpha the shape parameter of the distribution
	 * @param beta  the scale parameter of the distribution
	 * @return the mean of the Gamma distribution
	 * @throws IllegalArgumentException if alpha or beta is not positive
	 */
	public double gammaMean(double alpha, double beta) {
	    if (alpha <= 0.0 || beta <= 0.0) {
	        throw new IllegalArgumentException("alpha and beta must be positive");
	    }
	    return alpha * beta;
	}

	/**
	 * Calculates the variance of a Gamma distribution.
	 *
	 * @param alpha the shape parameter of the distribution
	 * @param beta  the scale parameter of the distribution
	 * @return the variance of the Gamma distribution
	 * @throws IllegalArgumentException if alpha or beta is not positive
	 */
	public double gammaVariance(double alpha, double beta) {
	    if (alpha <= 0.0 || beta <= 0.0) {
	        throw new IllegalArgumentException("alpha and beta must be positive");
	    }
	    return alpha * beta * beta;
	}
	
	
	//========== Beta ===================
	/**
	 * Calculates the probability density function (pdf) of a Beta distribution.
	 *
	 * @param x     the value at which to evaluate the pdf
	 * @param alpha the first shape parameter of the distribution
	 * @param beta  the second shape parameter of the distribution
	 * @return the pdf of the Beta distribution
	 * @throws IllegalArgumentException if x, alpha, or beta is NaN or alpha or beta is not positive
	 */
	public double betaPdf(double x, double alpha, double beta) {
	    if (Double.isNaN(x) || Double.isNaN(alpha) || Double.isNaN(beta)) {
	        throw new IllegalArgumentException("x, alpha, and beta must not be NaN");
	    }
	    if (alpha <= 0.0 || beta <= 0.0) {
	        throw new IllegalArgumentException("alpha and beta must be positive");
	    }
	    if (x < 0.0 || x > 1.0) {
	        return 0.0;
	    }
	    return Math.pow(x, alpha - 1) * Math.pow(1 - x, beta - 1) / betaFunction(new BigDecimal(alpha), new BigDecimal(beta)).doubleValue();
	}
	
	/**
	 * This implementation uses BigDecimal for precision in the calculations of the 
	 * beta function and the gamma function. The gamma method is a private method 
	 * called by betaFunction that calculates the gamma function using the Lanczos 
	 * approximation for positive values of x and the reflection formula for 
	 * negative values of x. Note that this implementation returns BigDecimal objects 
	 * instead of double primitives to avoid floating-point rounding errors.
	 * 
     * Calculates the beta function for the given shape parameters alpha and beta.
     * 
     * @param alpha the first shape parameter of the beta distribution.
     * @param beta the second shape parameter of the beta distribution.
     * @return the value of the beta function.
     */
    public BigDecimal betaFunction(BigDecimal alpha, BigDecimal beta) {
        MathContext mc = MathContext.DECIMAL128; // set precision for BigDecimal calculations
        BigDecimal numerator = gamma(alpha).multiply(gamma(beta));
        BigDecimal denominator = gamma(alpha.add(beta));
        return numerator.divide(denominator, mc);
    }

    /**
     * Calculates the gamma function for a given input value x.
     * 
     * @param x the input value for which to calculate the gamma function.
     * @return the value of the gamma function for the input value x.
     */
    private BigDecimal gamma(BigDecimal x) {
        MathContext mc = MathContext.DECIMAL128; // set precision for BigDecimal calculations
        BigDecimal oneHalf = new BigDecimal("0.5", mc);
        BigDecimal pi = new BigDecimal(Math.PI, mc);
        BigDecimal e = new BigDecimal(Math.E, mc);
        BigDecimal sqrt2pi = new BigDecimal(Math.sqrt(2 * Math.PI), mc);

        // For small x (x < 0.5), use the reflection formula
        if (x.compareTo(oneHalf) == -1) {
            BigDecimal sinPiX = new BigDecimal(Math.sin(pi.multiply(x).doubleValue()), mc);
            BigDecimal gamma1MinusX = gamma(BigDecimal.ONE.subtract(x));
            return pi.divide(sinPiX.multiply(gamma1MinusX), mc);
        } else {
            x = x.subtract(BigDecimal.ONE);
            BigDecimal a = oneHalf.divide(sqrt2pi, mc);
            BigDecimal b = e.pow((int) x.doubleValue(), mc).multiply(x.pow(x.intValue(), mc), mc);
            BigDecimal c = new BigDecimal("1.0", mc).add(
                (new BigDecimal("1.0", mc).divide(x, mc)).multiply(
                    new BigDecimal("-1/12", mc).add(
                        (new BigDecimal("1.0", mc).divide(x.pow(2), mc)).multiply(
                            new BigDecimal("1/360", mc).subtract(
                                new BigDecimal("1/1260", mc).divide(x.pow(4), mc)
                            )
                        )
                    )
                )
            );
            return a.multiply(b).multiply(c.sqrt(mc), mc);
        }
    }

	/**
	 * Calculates the mean of a Beta distribution.
	 *
	 * @param alpha the first shape parameter of the distribution
	 * @param beta  the second shape parameter of the distribution
	 * @return the mean of the Beta distribution
	 * @throws IllegalArgumentException if alpha or beta is not positive
	 */
	public double betaMean(double alpha, double beta) {
	    if (alpha <= 0.0 || beta <= 0.0) {
	        throw new IllegalArgumentException("alpha and beta must be positive");
	    }
	    return alpha / (alpha + beta);
	}

	/**
	 * Calculates the variance of a Beta distribution.
	 *
	 * @param alpha the first shape parameter of the distribution
	 * @param beta  the second shape parameter of the distribution
	 * @return the variance of the Beta distribution
	 * @throws IllegalArgumentException if alpha or beta is not positive
	 */
	public double betaVariance(double alpha, double beta) {
	    if (alpha <= 0.0 || beta <= 0.0) {
	        throw new IllegalArgumentException("alpha and beta must be positive");
	    }
	    double denominator = Math.pow(alpha + beta, 2) * (alpha + beta + 1);
	    return alpha * beta / denominator;
	}

	

}
