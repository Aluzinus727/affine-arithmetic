from solvers.af1 import AF1
from solvers.af2 import AF2

"""
Ejemplo del paper Extensions of Affine Arithmetic: Application to Unconstrained Global Optimization.
Página 12 del doc, seccion 3.3
"""
def run_paper_test():
    print('Paper test AF1')
    print('Resultado esperado: 25 + 0e_1 + e_2. Intervalo = [24, 26]')

    x = AF1(5, {1: 1})
    res = x * (10 - x)
    print(res, res.get_range())

    print('Paper test AF2')
    print('Resultado esperado: 25 + 0e_1 + 0e_2 + 0e_3 + 1e_4')

    x = AF2(5, {1: 1})
    res = x * (10 - x)
    print(res, res.get_range())

    res = 10*x - x**2
    print('res x**2', res, res.get_range())