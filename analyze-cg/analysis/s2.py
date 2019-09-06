import numpy as np

def calc_q(directors):
    """ Calculate the Q tensor from a list of directors

        [Qxx, Qxy, Qxz]
    Q = [Qyx, Qyy, Qyz]
        [Qzx, Qzy, Qzz]

    Parameters:
    -----------
    directors : list
        list of directors for each molecule

    Returns:
    --------
    Q : list
        3x3 array for Q tensor
    """

    directors = np.array(directors)
    Q = np.zeros((3, 3))
    Q = np.array([[np.sum(directors[:,i]*directors[:,j]*3) 
                    for j in range(3)] 
                    for i in range(3)])
    diag = np.array([len(directors)]*3)
    diag = np.diag(diag)
    Q = Q - diag
    Q /= 2 * len(directors)
    return Q

def calc_s2(q):
    """ Calculate the nematic order parameter (S2)
    S2 = dominant eigenvalue of the Q tensor

    Parameters:
    -----------
    q : list
        3x3 array for Q tensor

    Returns:
    --------
    float
    """

    w, _ = np.linalg.eig(q)
    return max(w)
