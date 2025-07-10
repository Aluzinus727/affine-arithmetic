import math
from solvers.base import Operable
from pyibex import Interval

class AF2(Operable):
    _epsilon = 1

    def __init__(
        self,
        base: float, 
        epsilons: dict={}, 
        non_affine: float=0, # x_(n+1)
        non_affine_p: float=0, #non affine positive, x_(n+2)
        non_affine_n: float=0 #non affine negative, x_(n+3)
    ):
        self.base = base
        self.epsilons = epsilons
        self.non_affine = non_affine
        self.non_affine_p = non_affine_p
        self.non_affine_n = non_affine_n

    def __str__(self):
        terms = [
            f"{'+' if coef >= 0 else '-'} {abs(coef):.1f}e_{idx}"
            for idx, coef in sorted(self.epsilons.items())  # 👈 orden por idx
        ]
        repr = f"{self.base:.1f}" + "".join(f" {term}" for term in terms) + f" + {self.non_affine}x_(n+1)"
        repr += f" + {self.non_affine_p}x_(n+2)"
        repr += f" + {self.non_affine_n}x_(n+3)"
        return repr

    @classmethod
    def _get_epsilon(cls):
        epsilon = cls._epsilon
        cls._epsilon += 1
        return epsilon

    @classmethod 
    def reset(cls):
        cls._epsilon = 1

    @classmethod
    def to_interval(cls, interval):
        a, b = interval
        base = (a + b) / 2
        error = (b - a) / 2

        new_epsilons = {
            AF2._get_epsilon(): error
        }

        return AF2(     
            base=base,
            epsilons=new_epsilons
        )

    @classmethod
    def _to_af2(cls, interval):
        a, b = interval
        base = (a + b) / 2
        error = (b - a) / 2
        non_affine = error
        new_epsilons = {}

        return AF2(     
            base=base,
            epsilons=new_epsilons,
            non_affine=non_affine
        )

    @classmethod
    def to_intervals(cls, intervals):
        converted_intervals = []
        for interval in intervals:
            new_af = cls.to_interval(interval)
            converted_intervals.append(new_af)

        return converted_intervals

    def get_range(self):
        total_deviation = sum(abs(coef) for coef in self.epsilons.values()) + self.non_affine
        #print(self.base, total_deviation, self.non_affine_n, self.non_affine_p)
        lb = self.base - total_deviation - self.non_affine_n # x_(n+2) debería tomar el valor 0 y x_(n+3) el valor -1
        ub = self.base + total_deviation + self.non_affine_p # x_(n+2) debería tomar el valor 1 y x_(n+3) el valor 0

        return (lb, ub)

    def __neg__(self):
        return self.__rsub__(0)

    def __add__(self, other):
        if not isinstance(other, (AF2, int, float)):
            return NotImplemented
        
        if isinstance(other, AF2):
            new_base = self.base + other.base
            new_epsilons = self.epsilons.copy()
            new_non_affine = self.non_affine + other.non_affine
            new_non_affine_p = self.non_affine_p + other.non_affine_p
            new_non_affine_n = self.non_affine_n + other.non_affine_n

            for idx, coef in other.epsilons.items():
                new_epsilons[idx] = new_epsilons.get(idx, 0.0) + coef

            return AF2(
                base=new_base,
                epsilons=new_epsilons,
                non_affine=new_non_affine,
                non_affine_p=new_non_affine_p,
                non_affine_n=new_non_affine_n
            )
        
        elif isinstance(other, (int, float)):
            self.base = self.base + other
            return self

    def __radd__(self):
        return self.__add__(other)

    def __sub__(self, other):
        if not isinstance(other, (AF2, int, float)):
            return NotImplemented

        if isinstance(other, AF2):
            new_base = self.base - other.base
            new_epsilons = self.epsilons.copy()
            new_non_affine = self.non_affine + other.non_affine

            for idx, coef in other.epsilons.items():
                new_epsilons[idx] = new_epsilons.get(idx, 0.0) - coef

            new_non_affine_p = self.non_affine_p + other.non_affine_n
            new_non_affine_n = self.non_affine_n + other.non_affine_p

        elif isinstance(other, (int, float)):
            new_base = self.base - other
            new_epsilons = self.epsilons.copy()
            new_non_affine = abs(self.non_affine)
            new_non_affine_p = self.non_affine_p
            new_non_affine_n = self.non_affine_n

        return AF2(
            base=new_base,
            epsilons=new_epsilons,
            non_affine=new_non_affine,
            non_affine_p=new_non_affine_p,
            non_affine_n=new_non_affine_n
        )

    def __rsub__(self, other):
        if not isinstance(other, (AF2, int, float)):
            return NotImplemented

        if isinstance(other, AF2):
            return other.__sub__(self)

        elif isinstance(other, (int, float)):
            new_base = other - self.base
            new_epsilons = self.epsilons.copy()
            for idx, coef in new_epsilons.items():
                new_epsilons[idx] = -coef

            new_non_affine = abs(self.non_affine)
            new_non_affine_p = self.non_affine_p
            new_non_affine_n = self.non_affine_n

            return AF2(
                base=new_base,
                epsilons=new_epsilons,
                non_affine=new_non_affine,
                non_affine_p=new_non_affine_p,
                non_affine_n=new_non_affine_n
            )

    def _k1(self, other):
        
        init = abs(self.base) * other.non_affine + abs(other.base) * self.non_affine
        sumatoria = 0
        for x_idx, x_coef in self.epsilons.items():
            for y_idx, y_coef in other.epsilons.items():
                if x_idx != y_idx:
                    sumatoria += abs(x_coef * y_coef)

            sumatoria += abs(x_coef * other.non_affine)
            sumatoria += abs(x_coef * other.non_affine_p)
            sumatoria += abs(x_coef * other.non_affine_n)
        
        #print(sumatoria)

        for y_idx, y_coef in other.epsilons.items():
            sumatoria += abs(y_coef * self.non_affine)
            sumatoria += abs(y_coef * self.non_affine_p)
            sumatoria += abs(y_coef * self.non_affine_n)
        
        #print(sumatoria)

        sumatoria += abs(self.non_affine * other.non_affine_p) + abs(self.non_affine * other.non_affine_n)
        sumatoria += abs(self.non_affine_p * other.non_affine) + abs(self.non_affine_p * other.non_affine_n)
        sumatoria += abs(self.non_affine_n * other.non_affine) + abs(self.non_affine_n * other.non_affine_p)
        
        #print(sumatoria)
        return init + sumatoria

    def _k2(self, other):
        x_0 = self.base
        y_0 = other.base

        if x_0 > 0 and y_0 > 0:
            init = x_0 * other.non_affine_p + y_0 * self.non_affine_p
        elif x_0 > 0 and y_0 < 0:
            init = x_0 * other.non_affine_p - y_0 * self.non_affine_n
        elif x_0 < 0 and y_0 > 0:
            init = -x_0 * other.non_affine_n + y_0 * self.non_affine_p
        elif x_0 < 0 and y_0 > 0:
            init = -x_0 * other.non_affine_n + y_0 * self.non_affine_n
        else:
            init = 0

        sumatoria = 0
        for x_idx, x_coef in self.epsilons.items():
            x_i = self.epsilons.get(x_idx, 0)
            y_i = other.epsilons.get(x_idx, 0)
            if x_i*y_i > 0:
                sumatoria += x_i * y_i

        if self.non_affine*other.non_affine > 0:
            sumatoria += self.non_affine * other.non_affine

        if self.non_affine_p*other.non_affine_p > 0:
            sumatoria += self.non_affine_p * other.non_affine_p

        if self.non_affine_n*other.non_affine_n > 0:
            sumatoria += self.non_affine_n * other.non_affine_n

        return init + sumatoria

    def _k3(self, other):
        x_0 = self.base
        y_0 = other.base

        if x_0 > 0 and y_0 > 0:
            init = x_0 * other.non_affine_n + y_0 * self.non_affine_n
        elif x_0 > 0 and y_0 < 0:
            init = x_0 * other.non_affine_n - y_0 * self.non_affine_p
        elif x_0 < 0 and y_0 > 0:
            init = -x_0 * other.non_affine_p + y_0 * self.non_affine_n
        elif x_0 < 0 and y_0 < 0:
            init = -x_0 * other.non_affine_p - y_0 * self.non_affine_p
        else:
            init = 0

        sumatoria = 0
        for x_idx, x_coef in self.epsilons.items():
            x_i = self.epsilons.get(x_idx, 0)
            y_i = other.epsilons.get(x_idx, 0)
            if x_i * y_i < 0:
                sumatoria += abs(x_i * y_i)

        if self.non_affine*other.non_affine < 0:
            sumatoria += abs(self.non_affine * other.non_affine)

        if self.non_affine_p*other.non_affine_p < 0:
            sumatoria += abs(self.non_affine_p * other.non_affine_p)

        if self.non_affine_n*other.non_affine_n < 0:
            sumatoria += abs(self.non_affine_n * other.non_affine_n)

        return init + sumatoria

    def __mul__(self, other):
        if not isinstance(other, (AF2, int, float)):
            return NotImplemented
    
        if isinstance(other, AF2):
            new_base = self.base * other.base

            # obtener todos los índices únicos
            all_indices = set(self.epsilons.keys()) | set(other.epsilons.keys())
            new_epsilons = {}
            for idx in all_indices:
                x_i = self.epsilons.get(idx, 0)
                y_i = other.epsilons.get(idx, 0)
                coef = self.base * y_i + x_i * other.base
                new_epsilons[idx] = coef

            new_non_affine = self._k1(other)
            new_non_affine_p = self._k2(other)
            new_non_affine_n = self._k3(other)

        elif isinstance(other, (int, float)):
            new_base = other * self.base
            new_epsilons = self.epsilons.copy()
            for idx, coef in self.epsilons.items():
                new_epsilons[idx] = self.epsilons[idx] * other

            if other > 0:                    
                new_non_affine = other * self.non_affine
                new_non_affine_p = other * self.non_affine_p
                new_non_affine_n = other * self.non_affine_n
            else:
                new_non_affine = abs(other) * self.non_affine
                new_non_affine_p = abs(other) * self.non_affine_p
                new_non_affine_n = abs(other) * self.non_affine_n

        return AF2(
            base=new_base,
            epsilons=new_epsilons,
            non_affine=new_non_affine,
            non_affine_p=new_non_affine_p,
            non_affine_n=new_non_affine_n
        )

    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                raise ValueError("División por cero")
            return self * (1 / other)

        elif isinstance(other, AF2):
            x_interval = self.get_range()
            y_interval = other.get_range()

            new_x_interval = Interval(min(x_interval[0], x_interval[1]), max(x_interval[0], x_interval[1]))
            new_y_interval = Interval[min(y_interval[0], y_interval[1]), max(y_interval[0], y_interval[1])]

            try:
                res = new_x_interval * (1 / new_y_interval)
                return self._to_af2((res.lb(), res.ub()))
            except:
                raise ValueError('Error __truediv__ AF2', x_interval, y_interval)

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            interval = self.get_range()
            x_interval = Interval(min(interval[0], interval[1]), max(interval[0], interval[1]))
            res = other / x_interval
            return self._to_af2((res.lb(), res.ub()))
            
        elif isinstance(other, AF2):
            x_interval = other.get_range()
            y_interval = self.get_range()

            new_x_interval = Interval(x_interval[0], x_interval[1])
            new_y_interval = Interval(y_interval[0], y_interval[1])

            try:
                res = new_x_interval * (1 / new_y_interval)
                return self._to_af2((res.lb(), res.ub()))
            except:
                raise ValueError('Error __rtruediv__ AF2', x_interval, y_interval)

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

        print('alpha, error, error_maximo', alpha, error, self)

        mult = alpha * self + error
        print('mult', mult, alpha*self)

        mult.non_affine += error_maximo
        return mult

    def __pow__(self, exponent):
        """
            #TODO: Reimplementar. No está ocupando e_(n+3) como corresponde.
            Según el paper, 10x - x^2 con x = 5 + e_1 debería dar 25 + x_(n+3) = 25 + [-1, 0] = (24, 25)
            Esta implementación da 24.5 + 0.5e_(n+1) = 24.5 + [-0.5, 0.5] = (24, 25)
            El intervalo está bien, pero la expresión no.
        """
        lb, ub = self.get_range()

        if exponent == 1:
            return self
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
        """
            Implementación según paper
            1. Se convierte x a intervalo
            2. Se calcula operación con aritmética intervalar tradicional
            3. Se vuelve a convertir a AF2.
        """
        lb, ub = self.get_range()

        if ub - lb >= 2 * math.pi:
            return self._to_af2((-1, 1))

        # Candidatos: extremos y puntos críticos donde cos(x) = ±1
        points = [lb, ub]

        n_start = math.floor(lb / math.pi) 
        n_end = math.floor(ub / math.pi)

        for n in range(n_start, n_end + 1):
            x = n * math.pi
            if lb <= x <= ub:
                points.append(x)

        values = [math.cos(x) for x in points]
        return self._to_af2((min(values), max(values)))

    def sin(self):
        print('Implementar sin. AF2')

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

        return self._to_af2((math.sqrt(lb), math.sqrt(ub)))

    def exp(self):
        """
            Implementación según paper
            1. Se convierte x a intervalo
            2. Se calcula operación con aritmética intervalar tradicional
            3. Se vuelve a convertir a AF2.
        """
        lb, ub = self.get_range()
        return self._to_af2((math.exp(lb), math.exp(ub)))
