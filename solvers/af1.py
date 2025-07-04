import math
from solvers.base import Operable
from pyibex import Interval

class AF1(Operable):
    _epsilon = 1

    def __init__(
        self, 
        base: float, 
        epsilons: dict={}, 
        non_affine: float=0
    ):
        self.base = base # x_0
        self.epsilons = epsilons # variables independientes
        self.non_affine = non_affine # x_(n + 1)

    def __str__(self):
        terms = [
            f"{'+' if coef >= 0 else '-'} {abs(coef):.1f}e_{idx}"
            for idx, coef in sorted(self.epsilons.items())  # 👈 orden por idx
        ]
        return f"{self.base:.1f}" + "".join(f" {term}" for term in terms) + f" + {self.non_affine}x_(n+1)"

    @classmethod
    def _get_epsilon(cls):
        epsilon = cls._epsilon
        cls._epsilon += 1
        return epsilon

    @classmethod 
    def reset(cls):
        cls._epsilon = 1

    @classmethod
    def _to_interval(cls, interval):
        a, b = interval
        base = (a + b) / 2
        error = (b - a) / 2

        new_epsilons = {
            AF1._get_epsilon(): error
        }

        return AF1(
            base=base,
            epsilons=new_epsilons
        )

    @classmethod
    def _to_af1(cls, interval):
        a, b = interval
        base = (a + b) / 2
        error = (b - a) / 2

        new_epsilons = {}
        non_affine = error

        return AF1(
            base=base,
            epsilons=new_epsilons,
            non_affine=non_affine
        )

    @classmethod
    def to_intervals(cls, intervals):
        converted_intervals = []
        for interval in intervals:
            new_af = cls._to_af1(interval)
            converted_intervals.append(new_af)

        return converted_intervals

    def get_range(self):
        total_deviation = sum(coef for coef in self.epsilons.values()) + self.non_affine
        return (self.base - total_deviation, self.base + total_deviation)

    def add_non_affine(self, new):
        self.non_affine += new
        self.non_affine = abs(self.non_affine)

    def __neg__(self):
        return self.__rsub__(0)

    def __add__(self, other):
        if not isinstance(other, (AF1, int, float)):
            return NotImplemented

        if isinstance(other, AF1):
            new_base = self.base + other.base
            new_epsilons = self.epsilons.copy()
            new_non_affine = self.non_affine + other.non_affine

            for idx, coef in other.epsilons.items():
                new_epsilons[idx] = new_epsilons.get(idx, 0.0) + coef

        elif isinstance(other, (int, float)):
            new_base = self.base + other
            new_epsilons = self.epsilons.copy()
            new_non_affine = abs(self.non_affine)

        return AF1(new_base, new_epsilons, new_non_affine)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if not isinstance(other, (AF1, int, float)):
            return NotImplemented

        if isinstance(other, AF1):
            new_base = self.base - other.base
            new_epsilons = self.epsilons.copy()
            new_non_affine = self.non_affine + other.non_affine

            for idx, coef in other.epsilons.items():
                new_epsilons[idx] = new_epsilons.get(idx, 0.0) - coef
        
        elif isinstance(other, (int, float)):
            new_base = self.base - other
            new_epsilons = self.epsilons.copy()
            new_non_affine = abs(self.non_affine)

        return AF1(new_base, new_epsilons, new_non_affine)

    def __rsub__(self, other):
        if not isinstance(other, (AF1, int, float)):
            return NotImplemented

        if isinstance(other, AF1):
            return other.__sub__(self)

        elif isinstance(other, (int, float)):
            new_base = other - self.base
            new_epsilons = self.epsilons.copy()
            for idx, coef in new_epsilons.items():
                new_epsilons[idx] = -coef

            new_non_affine = abs(self.non_affine)

            return AF1(new_base, new_epsilons, new_non_affine)

    def __mul__(self, other):
        if not isinstance(other, (AF1, int, float)):
            return NotImplemented

        if isinstance(other, AF1):
            new_base = self.base * other.base
            new_epsilons = {}

            # obtener todos los índices únicos
            all_indices = set(self.epsilons.keys()) | set(other.epsilons.keys())

            for idx in all_indices:
                x_i = self.epsilons.get(idx, 0)
                y_i = other.epsilons.get(idx, 0)
                coef = self.base * y_i + x_i * other.base
                new_epsilons[idx] = coef

            a = abs(self.base) * other.non_affine
            b = abs(other.base) * self.non_affine
            x_sum = abs(self.non_affine)
            y_sum = abs(other.non_affine)

            for idx, coef in self.epsilons.items():
                x_sum += abs(coef)

            for idx, coef in other.epsilons.items():
                y_sum += abs(coef)
                 
            new_non_affine = a + b + (x_sum * y_sum)

        if isinstance(other, (int, float)):
            new_base = other * self.base
            new_epsilons = self.epsilons.copy()
            new_non_affine = abs(other) * self.non_affine

            for idx, coef in self.epsilons.items():
                new_epsilons[idx] = self.epsilons[idx] * other

        return AF1(new_base, new_epsilons, new_non_affine)

    def __rmul__(self, other):
        return self.__mul__(other)
     
    def _apply_nonaffine_operation(self, func, other):
        """
            1. Convierte a intervalo
            2. Se calcula con interval arithmetic
            3. Se transforma el intervalo resultante a su forma afin nuevamente.
        """
        x_interval = self.get_range()
        x_min, x_max = x_interval
        
        # Aplicar la función al intervalo
        if x_min <= 0 <= x_max:
            raise ValueError("División por intervalo que contiene cero")
        if x_min > 0:
            result_min = other / x_max
            result_max = other / x_min
        else:  # x_max < 0
            result_min = other / x_max
            result_max = other / x_min

        # Crear nueva variable simbólica con el resultado
        return AffineForm.from_interval(result_min, result_max)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                raise ValueError("División por cero")
            return self * (1 / other)

        elif isinstance(other, AF1):
            x_interval = self.get_range()
            y_interval = other.get_range()

            new_x_interval = Interval[x_interval[0], x_interval[1]]
            new_y_interval = Interval[y_interval[0], y_interval[1]]

            try:
                res = new_x_interval * (1 / new_y_interval)
                return self._to_af1((res.lb(), res.ub()))
            except:
                raise ValueError('Error __truediv__ AF1', x_interval, y_interval)

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            interval = self.get_range()
            print('Interval', interval)
            x_interval = Interval(interval[0], interval[1])
            res = other / x_interval
            return self._to_af1((res.lb(), res.ub()))
            
        elif isinstance(other, AF1):
            x_interval = other.get_range()
            y_interval = self.get_range()

            new_x_interval = Interval(x_interval[0], x_interval[1])
            new_y_interval = Interval(y_interval[0], y_interval[1])

            try:
                res = new_x_interval * (1 / new_y_interval)
                return self._to_af1((res.lb(), res.ub()))
            except:
                raise ValueError('Error __rtruediv__ AF1', x_interval, y_interval)

    def chebyshev_approximation(self, func, df_proj):
        # Paso 1: Obtener el rango [a, b] de la forma afín actual
        a, b = self.get_range()
        
        # Paso 2: Evaluar la función en los extremos
        f_a = func(a)
        f_b = func(b)
        
        # Paso 3: Calcular el coeficiente α (pendiente de la línea secante)
        alpha = (f_b - f_a) / (b - a)
        
        # Paso 4: Encontrar u 
        u = df_proj(alpha)
        r = lambda x: alpha*x - alpha*a + f_a

        f_u = func(u)
        r_u = r(u)
        
        error = ((f_u + r_u) / 2) - (alpha * u)
        error_maximo =  abs(f_u - r_u) / 2

        mult = alpha * self + error
        mult.non_affine += error_maximo
        return mult

    def __pow__(self, exponent):
        """
            #TODO: Revisar. 
            Probablemente no ocupa cheby para aproximar.
            Se suma el error obtenido a la variable non_affine.
        """
        lb, ub = self.get_range()

        if exponent == 1:
            return AF1(self.center, self.linear_terms.copy(), self.error_term)
        elif exponent == 0:
            return AF1(1)
        elif exponent == 2:
            return self.chebyshev_approximation(lambda x: x**exponent, lambda a: a/2)
        elif lb >= 0 and exponent == 3:
            return self.chebyshev_approximation(lambda x: x**exponent, lambda a: math.sqrt(a/3))
        elif exponent == 4:
            return self.chebyshev_approximation(lambda x: x**exponent, lambda a: (a / 4) ** (1/3))
        elif exponent == 6:
            return self.chebyshev_approximation(lambda x: x**exponent, lambda a: (a / 6) ** (1/5))
        elif exponent == 7:
            return self.chebyshev_approximation(lambda x: x**exponent, lambda a: (a / 7) ** (1/6))

    def cos(self):
        lb, ub = self.get_range()

        if ub - lb >= 2 * math.pi:
            return self._to_af1((-1, 1))

        # Candidatos: extremos y puntos críticos donde cos(x) = ±1
        points = [lb, ub]

        n_start = math.floor(lb / math.pi) 
        n_end = math.floor(ub / math.pi)

        for n in range(n_start, n_end + 1):
            x = n * math.pi
            if lb <= x <= ub:
                points.append(x)

        values = [math.cos(x) for x in points]
        return self._to_af1((min(values), max(values)))

    def sin(self):
        print('Implementar sin. AF1')

    def sqrt(self):
        """
            Implementación según paper
            1. Se convierte x a intervalo
            2. Se calcula operación con aritmética intervalar tradicional
            3. Se vuelve a convertir a AF2.
        """
        lb, ub = self.get_range()
        if lb < 0:
            raise ValueError('Raiz de valor negativo', lb)

        return self._to_af1((math.sqrt(lb), math.sqrt(ub)))
    
    def exp(self):
        """
            Implementación según paper
            1. Se convierte x a intervalo
            2. Se calcula operación con aritmética intervalar tradicional
            3. Se vuelve a convertir a AF2.
        """
        lb, ub = self.get_range()
        return self._to_af1((math.exp(lb), math.exp(ub)))