import numpy as np
from scipy.signal import gaussian, convolve


class GaussianSmoothing:
    def __init__(self, raw, sigma=3, length=10, axis='row', decimal_out=2):
        """
        :param raw: raw signal (1d, e.g., spike trains or binned spike counts)
        :param sigma: standard deviation of the gaussian kernel
        :param length: width of the gaussian kernel
        """
        assert raw.ndim <= 2  # designed to deal with 2-d arrays
        self.raw = np.array(raw)
        self.sigma = sigma
        self.length = length
        self.axis = axis
        self.decimal_out = decimal_out
        # self.norm_factor = 1/(self.sigma * np.sqrt(2 * np.pi))

    def get_gaussian_kernel(self):
        """
        Uses scipy.signal.gaussian to generate a gaussian kernel
        :return:
        """
        g = gaussian(self.length, self.sigma)
        g_norm = g/np.sum(g)  # normalize
        return g_norm

    def conv(self):
        out = np.empty(self.raw.shape)
        out[:] = np.nan
        kernel = self.get_gaussian_kernel()
        if self.axis == 'row':
            for i in range(self.raw.shape[0]):
                out[i, :] = convolve(self.raw[i, :], kernel, mode='same')
        elif self.axis == 'column':
            for i in range(self.raw.shape[1]):
                out[:, i] = convolve(self.raw[:, i], kernel, mode='same')
        return np.around(out, self.decimal_out)
