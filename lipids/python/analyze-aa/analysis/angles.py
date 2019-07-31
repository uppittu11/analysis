import numpy as np
### This module is used to calculate things with angles


def calc_angle(vec1, vec2):
    """ Calculate the angles between two lists of unit vectors
    theta = arccos( dot(vec1*vec2) )

    Parameters:
    -----------
    vec1 : list
        list of unit vectors
    vec2 : list
        list of unit vectors

    Returns:
    --------
    theta : list
    """

    dot = np.sum(np.array(vec1)*np.array(vec2), axis=1)
    theta = np.arccos(dot)*180/np.pi
    mask = (theta > 90).astype(int)
    theta = (180*mask) - theta * (-1*mask)
    return theta

def calc_direction_vector(coord1, coord2):
    """ Calculate the vector between two points
    vec = (coord1 - coord2) / norm(coord1-coord2)

    Parameters:
    -----------
    coord1 : list
        list of unit coordinates
    coord2 : list
        list of unit coordinates

    Returns:
    --------
    vector : list
    """

    vector = np.array(coord2)-np.array(coord1)
    magnitude = np.linalg.norm(axis=1)
    vector = vector/magnitude[:,np.newaxis]
    return vector
