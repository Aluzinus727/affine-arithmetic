from pyibex import Interval
from solvers.paper import AFPaper
from solvers.af2 import AF2

"""
Ejemplo propio del paper que se está escribiendo de AF revisada.
"""
interval_1 = Interval(0, 2)
interval_2 = Interval(0, 2)

x = AFPaper(a=[1, 2], x=[interval_1, interval_2], error=Interval(-1, 1))
y = AFPaper(a=[2, 3], x=[interval_1, interval_2], error=Interval(0, 2))

print('Ejemplo con AF revisada')
print(f"Intervalo 1: {interval_1}. X: {x}")
print(f"Intervalo 2: {interval_2}. Y: {y}")

res_afr = x*y
interval_afr = res_afr.to_interval()
rango_afr = interval_afr[1] - interval_afr[0]
print('Resultado ', res_afr, interval_afr, rango_afr, end='\n\n')

x = AF2.to_interval([-1, 7])
y = AF2.to_interval([0, 12])

res = x*y
interval = res.get_range()
rango = interval[1] - interval[0]
print('Ejemplo con AF2')
print('Resultado ', res, interval, rango)
