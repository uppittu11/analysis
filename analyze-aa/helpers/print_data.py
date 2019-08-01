import numpy as np
import scipy.stats as stats
import pickle

# Load in the results file
with open('results.p', 'rb') as f:
    results = pickle.load(f)

print(results)

# Check if the results are for a bilayer (top and bottom) or
# for a multilayer (no distinction)
if len(results[0]) <= 5:
    tilt = results[:,:2]
    print("Tilt :   {:.5} +/- {:.5}".format(np.mean(tilt[:,0]), np.linalg.norm(tilt[:,1])))
    s2 = results[:,2]
    print("S2 :     {:.5} +/- {:.5}".format(np.mean(s2), stats.sem(s2)))
    apl = results[:,3]
    print("APL :    {:.5} +/- {:.5}".format(np.mean(apl), stats.sem(apl)))
    height = results[:,4]
    print("Height : {:.5} +/- {:.5}".format(np.mean(height), stats.sem(height)))
