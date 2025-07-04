"""
Este archivo contiene la clase que guardará la memoria de las cajas.
"""
class Affine:
    x_ref = None
    envs = None

    @classmethod
    def reset(cls):
        cls.x_ref = None
        cls.envs = dict()
