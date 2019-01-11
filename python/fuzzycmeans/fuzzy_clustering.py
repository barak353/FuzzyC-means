import numpy as np
import random


Epsilon = 0.00001


class FCM:
    """
    Pseudo code:
        1) initialized W matrix
        2) compare the centers of the clusters and if it not different from the old one do:
            3) update W matrix
            4) compute the centers of the clusters
    """

    def __init__(self, numOfClusters=2, p=9):

        self.numOfClusters = numOfClusters
        # The fuzziness.
        self.p = p
        # W is the weights matrix.
        self.W = None
        # Centroids of each cluster.
        self.clustersCenters = None

    def InitMembershipMatrix(self, numOfPoints):
        """
        This function initialize the weights matrix.
        :param numOfPoints: number of points define the size of the W matrix.
        """
        self.W = np.zeros((numOfPoints, self.numOfClusters))
        for i in range(0, numOfPoints):
            rowSum = 0.0
            for j in range(0, self.numOfClusters - 1):
                randNum = random.random()
                randNum = round(randNum, 2)
                # The sum of the weights can't be larger than 1.
                if randNum + rowSum <= 1.0:
                    self.W[i][j] = randNum
                    rowSum += self.W[i][j]
            # If it's the last iteration
            self.W[i][self.numOfClusters-1] = 1.0 - rowSum

    def ComputeCentroids(self, X):
        """
        This function calculate the centroids of each cluster.
        :param X: set of points
        :return: List of centers.
        """
        n = X.shape[0]
        m = X.shape[1]
        centers = []
        # For each cluster
        for c in range(0, self.numOfClusters):
            # Weight
            WeightsOfXvec = np.zeros(m)
            sumOfWeights = 0.0
            for i in range(0, n):
                # (numerator\denominator) = sigma[0-n]((W[i][c])^2) * X[i]\ sigma[0-n]((W[i][c])^2)
                denominator = (self.W[i][c] ** self.p)
                numerator = denominator * X[i]
                WeightsOfXvec += numerator
                sumOfWeights += denominator
            if sumOfWeights == 0:
                # To prevent division in zero.
                sumOfWeights = 0.000001
            centers.append(WeightsOfXvec / sumOfWeights)
        # Update list.
        self.clustersCenters = centers
        return centers

    def EuclideanDistance(self, x, c):
        """
        Compute the Euclidean distance between X and C
        :param x: set of X points
        :param c: the center or a cluster
        :return: the distance
        """
        dis = 0.0
        for i in range(0, len(x)):
            dis += (x[i]-c[i]) ** 2
        return dis

    def UpdateWeights(self, X):
        """
        update the weights matrix
        :param X: data points
        """
        for i in range(0, X.shape[0]):
            for c in range(0, len(self.clustersCenters)):
                self.W[i][c] = self.ComputeWeight(X, i, c)

    def ComputeWeight(self, X, i, c):
        """
        :param X: data points
        :param i: specific element in x
        :param c: specific cluster
        :return: weight of the given element
        """
        numerator = self.EuclideanDistance(X[i], self.clustersCenters[c])
        # To prevent from return infinity.
        if numerator == 0:
            numerator = 1.0 - Epsilon
        numerator = (1 / numerator) ** (1.0 / (self.p - 1))
        sum = 0.0
        # (numerator\denominator) = (1 / dist(xi,cj)^(1\p-1)) \ sigma[0-numOfClusters](1/dist(xi,cj)^(1\p-1))
        for cluster in self.clustersCenters:
            denominator = self.EuclideanDistance(cluster, X[i])
            if denominator == 0.0:
                denominator = 1 / Epsilon
            sum += ((1 / denominator) ** (1.0 / (self.p - 1)))
        # To prevent from return infinity.
        if sum == 0:
            return 1.0 - Epsilon
        return (numerator / sum)

    def Algorithm(self, X):
        """
        :param X: set of points
        """
        computeCenters = False
        if self.clustersCenters is None:
            computeCenters = True
        if self.W is None:
            numOfPoints = X.shape[0]
            self.InitMembershipMatrix(numOfPoints)
        centersList = []
        i = 0
        while computeCenters and i < 50:
            i += 1
            centers = self.ComputeCentroids(X)
            centersList.append(centers)
            centersList = np.array(centersList).tolist()
            self.UpdateWeights(X)
            if centersList.__len__() >= 3:
                if any(x in centersList[centersList.__len__() - 1] for x in centersList[centersList.__len__() - 3]):
                    computeCenters = False
        return self.W