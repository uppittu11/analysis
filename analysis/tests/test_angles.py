import pytest
import numpy as np
import analysis
from analysis.tests.basetest import BaseTest

class TestAngles(BaseTest):
    @pytest.fixture
    def sample_coords(self):
        coord1 = [[1, 0, 0],
                  [0, -1, 0],
                  [1, 1, 1],
                  [-1, 2, 4]]
        coord2 = [[0, 1, 0],
                  [0, 1, 0],
                  [1, 2, 1],
                  [2, 3, -1]]
        coord1 = np.array(coord1, dtype=float)
        coord2 = np.array(coord2, dtype=float)

        return coord1, coord2
    
    @pytest.fixture
    def sample_vectors(self):
        vec1 = [[0, 0, 1],
                [-1, 0, 0],
                [1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]]
        vec2 = [[0, 0, 1],
                [1, 0, 0],
                [0, 1/np.sqrt(2), 1/np.sqrt(2)]]
        vec1 = np.array(vec1)
        vec2 = np.array(vec2)

        return vec1, vec2


    def test_calc_direction_vector(self, sample_coords):
        coord1, coord2 = sample_coords
        testfn = analysis.angles.calc_direction_vector
        result = testfn(coord1, coord2)
        expected = [[-1, 1, 0],
                    [0, 2, 0],
                    [0, 1, 0],
                    [3, 1, -5]]
        expected = np.array(expected, dtype=float)
        magnitudes = np.sqrt(expected[:,0] **2 +
                             expected[:,1] **2 +
                             expected[:,2] **2)
        expected /= magnitudes[:,np.newaxis]
        assert np.allclose(result, expected)

    def test_calc_angle(self, sample_vectors):
        vec1, vec2 = sample_vectors
        testfn = analysis.angles.calc_angle
        result = testfn(vec1, vec2)
        expected = [0.0,
                    0.0,
                    35.2643896827547]
        assert np.allclose(result, expected)
        