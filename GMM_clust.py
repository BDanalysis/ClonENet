import numpy as np
import random

def colors(n):
    ret = []
    for i in range(n):
        ret.append((random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)))
    return ret

class GMM(object):
    """Gaussian Mixture Model 
    """

    def __init__(self, data, K):
        """
        K: the number of gaussian models
        alpha: the weight for corresponding gaussian model
        mu: the vector of means
        sigma2: the vector of variances
        N: the number of samples
        K: the number of gaussian models
        """
        self.data = data
        self.K = K
        self.alpha = (np.ones(K) / K)
        self.mu = (np.arange(K) - K // 2) * (data.max() - data.min()) / K
        self.sigma2 = (np.ones(K))
        self.N = len(data)
        self.gamma = np.ones((self.N, K)) / K

    def phi(self):
        # phi.shape(K, N)
        phi = (1 / np.sqrt(2 * np.pi * self.sigma2.reshape(self.K, 1)) * np.exp(- (self.data - self.mu.reshape(self.K, 1)) ** 2 / (2 * self.sigma2.reshape(self.K, 1))))
        return phi

    def fit(self):
        sigma2_ = self.sigma2
        mu_ = self.mu
        k = 0
        while True:
            # gamma.shape(N, K)
            self.gamma = (0.1 * self.gamma + 0.9 * self.phi().T * self.alpha / (self.phi().T * self.alpha).sum(axis=1).reshape(self.N, 1))
            # mu.shape(1, K)
            self.mu = (0.1 * self.mu + 0.9 * np.matmul(self.data, self.gamma) / self.gamma.sum(axis=0))
            # sigma2.shape(1,K)
            self.sigma2 = (0.1 * self.sigma2 + 0.9 * (self.gamma * (self.data.reshape(self.N, 1) - self.mu) ** 2).sum(axis=0) / self.gamma.sum(axis=0))
            # alpha.shape(1, K)
            self.alpha = (0.1 * self.alpha + 0.9 * self.gamma.sum(axis=0) / self.N)
            if (np.sum((self.mu - mu_) ** 2) + np.abs(self.sigma2 - sigma2_).sum()) < 10 ** (-10):
                break
            mu_ = self.mu
            sigma2_ = self.sigma2
            k += 1
            #print(self.gamma.argmax(axis=1))

        return self.gamma.argmax(axis=1)