import pytest


class BaseTest:
    def get_fn(self, filename):
        import os.path

        full_path = os.path.join(__file__, "include", filename)
        return full_path
