## !/usr/bin/python
#     pylint: disable=E1101
from numpy.random import permutation, uniform, randint
from ea import DEcrossover
from DE import EAresult

from scipy.optimize import fmin_l_bfgs_b

import numpy as np
import time

import SHADE
from mts import mtsls


dict = {'mts': 0.9, 'grad': 0.1}
value = np.random.choice(list(dict.keys()))
print(value)
