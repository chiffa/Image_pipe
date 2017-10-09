import os
import unittest
os.environ['UNITTESTING'] = 'True'


class HooksConfigTest(unittest.TestCase):

    def test_hooks(self):
        self.assertTrue(True)

if __name__ == "__main__":
    my_list = []
    unittest.main()
