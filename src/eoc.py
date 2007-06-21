import pylinear.array as num
import pylinear.toybox as toybox




def estimate_order_of_convergence(abscissae, errors):
    """Assuming that abscissae and errors are connected by a law of the form

    error = constant * abscissa ^ (-order),

    this function finds, in a least-squares sense, the best approximation of
    constant and order for the given data set. It returns a tuple (constant, order).
    Both inputs must be PyLinear vectors.
    """

    assert len(abscissae) == len(errors)
    if len(abscissae) <= 1:
        raise RuntimeError, "Need more than one value to guess order of convergence."

    coefficients = toybox.fit_polynomial(num.log10(abscissae), num.log10(errors), 1)
    return 10**coefficients[0], -coefficients[1]


  




class EOCRecorder:
    def __init__(self):
        self.history = []

    def add_data_point(self, abscissa, error):
        self.history.append((abscissa, error))

    def estimate_order_of_convergence(self, gliding_mean = None):
        abscissae = num.array([ a for a,e in self.history ])
        errors = num.array([ e for a,e in self.history ])

        size = len(abscissae)
        if gliding_mean is None:
            gliding_mean = size

        data_points = size - gliding_mean + 1
        result = num.zeros((data_points, 2), num.Float)
        for i in range(data_points):
            result[i,0], result[i,1] = estimate_order_of_convergence(
                abscissae[i:i+gliding_mean], errors[i:i+gliding_mean])
        return result

    def write_gnuplot_file(self, filename):
        outfile = file(filename, "w")
        for absc, err in self.history:
            outfile.write("%f %f\n" % (absc, err))
        result = self.estimate_order_of_convergence()
        const = result[0,0]
        order = result[0,1]
        outfile.write("\n")
        for absc, err in self.history:
            outfile.write("%f %f\n" % (absc, const * absc**(-order)))

