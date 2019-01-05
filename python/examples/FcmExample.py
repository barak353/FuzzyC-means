import numpy as np
import sys
sys.path.append('..')
from fuzzycmeans import FCM
import time

def example():
    X = np.array([[1, 1, 1], [1, 2, 2], [2, 2, 3], [9, 10, 11], [10, 10, 10], [10, 9, 9], [9, 9, 9], [20, 20, 20]])
    fcm = FCM(numOfClusters=15)
    fcm.Algorithm(X)

start_time = time.time()
example()
print("--- %s seconds ---" % (time.time() - start_time))
