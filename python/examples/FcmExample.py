import numpy as np
import sys
sys.path.append('..')
from fuzzycmeans import FCM
import matplotlib.pyplot as plt
import time

def example():
    X = np.array([[1, 1], [1, 2], [2, 2], [9, 10], [10, 10], [10, 9], [9, 9], [20,20]])
    fcm = FCM(numOfClusters=8)
    return fcm.Algorithm(X)

start_time = time.time()
res = example()
print("--- %s seconds ---" % (time.time() - start_time))
print("result:")
print(res)
plt.matshow(np.asarray(res))
plt.colorbar()
plt.show()