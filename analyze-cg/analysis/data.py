import numpy as np
import pickle

def load_results(filename, convert_to_numpy=True):
    with open(filename, 'rb') as f:
        results = pickle.load(f)
    if convert_to_numpy:
        results = _to_dict(results)

    return results
    
def _to_dict(results, fields=None):
    """ Convert a results list of dicts to a dict of lists
    Results are outputted as a list of length n_frames containing 
    dicts with keys corresponding to each result type. This function
    converts that format to a dict of numpy arrays for each result type
    in which all frames are in one data structure

    Parameters:
    -----------
    results : list
        Results in a list of dicts format
    fields : list
        Which keys in the dict to convert and output

    Returns:
    --------
    results : dict
        results in a dict of np.ndarrays format
    """

    new_results = dict()

    if not fields:
        fields = results[0].keys
    
    for field in fields:
        temp = [result[field] for result in results]
        temp = np.array(temp)
        new_results.update({field : temp})
        
    return new_results

