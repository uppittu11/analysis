import numpy as np
import scipy.stats as stats
import pickle

# Load in the results file
with open('results.p', 'rb') as f:
    results = pickle.load(f)

# Check if the results are for a bilayer (top and bottom) or
# for a multilayer (no distinction)
if len(results[0]) <= 5:
    tilt = np.array([result[0] for result in results])
    print("Tilt :   {:.5f} +/- {:.5f}".format(np.mean(tilt), stats.sem(tilt, axis=None)))
    s2 = np.array([result[1] for result in results])
    print("S2 :     {:.5f} +/- {:.5f}".format(np.mean(s2), stats.sem(s2)))
    apl = np.array([result[2] for result in results]) * 36
    print("APL :    {:.5f} +/- {:.5f}".format(np.mean(apl), stats.sem(apl)))
    height = np.array([result[3] for result in results]) * 6
    print("Height : {:.5f} +/- {:.5f}".format(np.mean(height), stats.sem(height, axis=None)))
