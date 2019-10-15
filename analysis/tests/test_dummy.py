import pytest
import analysis
from analysis.tests.basetest import BaseTest

class TestDummy(BaseTest):
    def test_dummy(self):
        assert 10 == 11.0

