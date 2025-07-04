import math
import torch
from typing import List
from pyibex import Interval, IntervalVector
from memory import Affine
from utils import circular_queue

class IntervalTensor:
    def __init__(self, lower, upper=None):
        if upper is None: upper=lower

        if not isinstance(lower, torch.Tensor):
            lower = torch.as_tensor(lower)  # no seteas requires_grad explícitamente
        if not isinstance(upper, torch.Tensor):
            upper = torch.as_tensor(upper)
            
        if isinstance(lower, torch.Tensor):
            self.lower = lower
        else:
            self.lower = torch.tensor(lower, dtype=torch.float32, requires_grad=False)

        if isinstance(upper, torch.Tensor):
            self.upper = upper
        else:
            self.upper = torch.tensor(upper, dtype=torch.float32, requires_grad=False)

        if self.lower > self.upper:
            self.lower, self.upper = self.upper, self.lower

    def __str__(self):
        return f'IntervalTensor({self.lower}, {self.upper})'

    def __repr__(self):
        return self.__str__()

    def __truediv__(self, other):
        res = self * other**(-1)
        return res

    def __pow__(self, exponent):
        if isinstance(exponent, int):
          if exponent==2:
            if self.lower <=0 and self.upper >= 0:
               new_lower = torch.tensor(0.0, dtype=torch.float32, requires_grad=False)
               new_upper = torch.max(self.lower**2, self.upper**2)
            elif self.upper <0:
               new_lower = self.upper ** 2.
               new_upper = self.lower ** 2.
            else:
               new_lower = self.lower ** 2.
               new_upper = self.upper ** 2.

        else:
            raise ValueError(f"exponent={exponent} is not implemented")

        return IntervalTensor(new_lower, new_upper)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if isinstance(other, IntervalTensor):
            # Caso 1: lower < 0 y upper < 0
            if self.lower < 0 and self.upper < 0:
                new_lower = self.upper * other.upper
                new_upper = self.lower * other.lower

            # Caso 2: lower < 0 y upper > 0
            elif self.lower < 0 and self.upper > 0:
                new_lower = self.lower * other.upper
                new_upper = self.upper * other.upper

            # Caso 3: lower > 0 y upper > 0
            elif self.lower > 0 and self.upper > 0:
                new_lower = self.lower * other.lower
                new_upper = self.upper * other.upper

            return IntervalTensor(new_lower, new_upper)
        elif isinstance(other, (int, float, torch.Tensor)):
            # Escalar o tensor regular multiplicando un intervalo
            new_lower = self.lower * other
            new_upper = self.upper * other
            return IntervalTensor(new_lower, new_upper)
        elif isinstance(other, LinearExpression):
            other.__mul__(self)
        else:
            raise NotImplementedError("Unsupported multiplication type")

    def __add__(self, other):
        if type(other) == IntervalTensor:
            new_lower = self.lower + other.lower
            new_upper = self.upper + other.upper
            return IntervalTensor(new_lower, new_upper)
        elif isinstance(other, (int, float, torch.Tensor)):
            new_lower = self.lower + other
            new_upper = self.upper + other
            return IntervalTensor(new_lower, new_upper)
        else:
            raise NotImplementedError("Unsupported addition type")

    def __neg__(self):
        return IntervalTensor(-self.upper,-self.lower)

    def __sub__(self, other):
        return self + (-other)
    
    def __contains__(self, item):
        return torch.all((self.lower <= item) & (item <= self.upper)).item()
    
    def __iter__(self):
        yield self.lower
        yield self.upper

    def mid(self):
        return (self.lower+self.upper)/2.

    def diam(self):
        return self.upper - self.lower

    def ub(self):
        # Calcula el máximo teórico del intervalo
        return self.upper

    def lb(self):
        # Calcula el máximo teórico del intervalo
        return self.lower

    def __str__(self):
       return f'IntervalTensor([{self.lower.item()} {self.upper.item()}])'


#Definición de la Clase Expresión Lineal
class LinearExpression:
    def __init__(
        self, 
        a: List, 
        b: int, 
        x: List[IntervalTensor]
    ):
        # f(x)= a1*x1 + a2*x2 + ... + b
        self.a = a
        self.b = b

        self.x = x

    #Llamada como función
    def __call__(self, x=None):
        if x:
            return sum(ai * xi for ai, xi in zip(self.a, x)) + self.b
        else:
            return sum(ai * xi for ai, xi in zip(self.a, self.x)) + self.b

    @classmethod
    def reset(cls):
        Affine.reset()

    @classmethod
    def to_intervals(cls, intervals):
        vars = [Interval(min(inter[0], inter[1]), max(inter[0], inter[1])) for inter in intervals]
        new_x = IntervalVector(vars)
        tensor_vars = [IntervalTensor(min(inter[0], inter[1]), max(inter[0], inter[1])) for inter in intervals]
       

        Affine.x_ref = tensor_vars
        Affine.envs = dict()

        n = new_x.size()
        linear_vars = []
        for i in range(n):
            a = [1 if j == i else 0 for j in range(n)]
            linear_vars.append(LinearExpression(a,0, tensor_vars))

        envs = []
        for i, var in enumerate(linear_vars):
            Affine.envs[i] = LinearEnvelope(is_var=True, linear_expr=var)
            envs.append(Affine.envs[i])

        return linear_vars

    def get_range(self):
        pass

    def cheby_convex(self, f=lambda x: x**2, df_proj=lambda m: m/2, convex=False, concave=False):
        x = self()
        lb, ub = x.lb(), x.ub()

        requires_grad = lb.requires_grad or ub.requires_grad
        data = circular_queue.get_next()

        initial_m = (f(ub) - f(lb)) / x.diam()
        if data:
            alfa = data['alfa']
        else:
            alfa = torch.tensor(1, dtype=torch.float64, requires_grad=True)
            m_data = {
                'alfa': alfa,
            }
            circular_queue.add(m_data)

        m = alfa * initial_m

        # Identificar el Punto de Tangencia xp
        # Necesitamos resolver df(xp) = m
        xp = df_proj(m)

        if convex: 
            if m >= initial_m:
                upper = m*self
                upper.b += (f(lb) - m * lb)
            else:
                upper = m*self 
                upper.b += (f(ub) - m * ub)

            lower = m*self
            lower.b += (f(xp) - m * xp)
        elif concave:
            if m >= initial_m:
                lower = m*self 
                lower.b += (f(ub) - m * ub)
            else:
                lower = m*self 
                lower.b += (f(lb) - m * lb)

            upper = m*self
            upper.b += (f(xp) - m * xp)

        return LinearEnvelope(lower = lower, upper= upper)

    def __pow__(self, exponent):
        if exponent == 2:
            return self.cheby_convex(f=lambda x: x**2, df_proj=lambda m: m/2, convex=True)
        elif exponent == 3:
            return self.cheby_convex(f=lambda x: x**3, df_proj=lambda m: math.sqrt(m/3), convex=True)
        elif exponent == 4:
            return self.cheby_convex(f=lambda x: x**4, df_proj=lambda m: (m / 4) ** (1 / 3), convex=True)
        elif exponent == 5:
            return self.cheby_convex(f=lambda x: x**5, df_proj=lambda m: (m / 5) ** (1 / 4), convex=True)
        elif exponent == 6:
            return self.cheby_convex(f=lambda x: x**6, df_proj=lambda m: (m / 6) ** (1 / 5), convex=True)
        elif exponent == 7:
            return self.cheby_convex(f=lambda x: x**7, df_proj=lambda m: (m / 7) ** (1 / 6), convex=True)
        
    def sqrt(self):
        return self.cheby_convex(
            f=lambda x: math.sqrt(x),
            df_proj=lambda m: 1 / (4 * math.pow(m, 2)),
            concave=True)


    def exp(self):
        return self.cheby_convex(
            f=lambda x: math.exp(x),
            df_proj=lambda m: math.log(m),
            convex=True)

    def cos(self):
        x = self()
        lb, ub = x.lb(), x.ub()

        if (ub - lb) > math.pi:
            raise ValueError("El intervalo es demasiado amplio, no se puede determinar convexidad/concavidad.", lb, ub)

        def calculate_n_convex(lb, ub):
            n_lb = math.floor((lb - math.pi / 2) / (2 * math.pi))
            n_ub = math.floor((ub - math.pi / 2) / (2 * math.pi))

            if n_lb == n_ub:
                return n_lb
            raise ValueError("cos, intervalo pasa por múltiples regiones convexas", lb, ub)

        def calculate_n_concave(lb, ub):
            n_lb = math.floor((lb + math.pi / 2) / (2 * math.pi))
            n_ub = math.floor((ub + math.pi / 2) / (2 * math.pi))

            if n_lb == n_ub:
                return n_lb
            raise ValueError("cos, ntervalo pasa por múltiples regiones cóncavas", lb, ub)

        # Verificar si está en una región convexa
        if math.cos(lb) <= 0 and math.cos(ub) <= 0:
            n = calculate_n_convex(lb, ub)
            inter_a = math.pi / 2 + 2*n*math.pi
            inter_b = (3 * math.pi) / 2 + 2*n*math.pi
            return self.cheby_convex(
                f=lambda x: math.cos(x),
                df_proj=lambda m: math.pi - math.asin(-m) + 2 * n * math.pi,
                convex=True
            )

        # Verificar si está en una región cóncava
        elif math.cos(lb) >= 0 and math.cos(ub) >= 0:
            n = calculate_n_concave(lb, ub)
            inter_a = -math.pi / 2 + 2*n*math.pi
            inter_b = math.pi / 2 + 2*n*math.pi
            return self.cheby_convex(
                f=lambda x: math.cos(x),
                df_proj=lambda m: math.asin(-m) + 2*n*math.pi,
                concave=True
            )

        raise ValueError("cos not implemented", lb, ub)

    #Suma de una expresión lineal con otra expresión lineal
    def __add__ (self, other):
        if isinstance(other, (int, float, torch.Tensor)):
            new_a = self.a
            new_b = self.b + other
            return LinearExpression(new_a, new_b, self.x)
        elif isinstance(other, (IntervalTensor)):
            new_a = self.a
            new_b = self.b + other
            return LinearExpression(new_a, new_b, self.x)
        elif isinstance(other, LinearExpression):
            new_a = [ai + bi for ai, bi in zip(self.a, other.a)]
            new_b = self.b + other.b
            return LinearExpression(new_a, new_b, self.x)
        elif isinstance(other, LinearEnvelope):
            return LinearEnvelope(lower = other.lower + self, upper=other.upper + self)

    #Multiplicación
    def __mul__(self, other):
        if isinstance(other, (int, float, torch.Tensor)):
            new_a = [ai * other for ai in self.a]
            new_b = self.b * other
            return LinearExpression(new_a, new_b, self.x)

        if isinstance(other, LinearExpression):
            v = self; w = other
            if(len(v.a) != len(w.a)): raise NotImplementedError("Las expresiones lineales tienen diferentes cantidades de variables.")

            env = LinearEnvelope() # 0

            # Sumatoria de (a_0 * b_i + b_0 * a_i) x_i
            for i, (ai, bi) in enumerate(zip(v.a, w.a)):
                if ai != 0 or bi != 0:
                    env += (v.b * bi + w.b * ai) * Affine.envs[i]

            # Sumatoria de a_i * b_i * x_i^2
            for i, (ai, bi) in enumerate(zip(v.a, w.a)):
                if ai != 0 and bi != 0:
                    env += (ai * bi) * Affine.envs.setdefault((i,i), Affine.envs[i] ** 2)

            # Sumatoria de a_i * b_j * x_i * x_j
            for i, ai in enumerate(v.a):
                for j, bj in enumerate(w.a):
                    if i != j and ai != 0 and bj != 0:
                        env_mult = Affine.envs.setdefault((i, j), Affine.envs[i] * Affine.envs[j])
                        Affine.envs.setdefault((j, i), env_mult)
                        env += ai * bj * env_mult

            # Agregar a_0 * b_0
            data = circular_queue.get_next()
            if data:
                alfa = data['alfa']
            else:
                alfa = torch.tensor(1, dtype=torch.float64, requires_grad=True)
                m_data = {
                    'alfa': alfa,
                }
                circular_queue.add(m_data)
            new_env = env * alfa
            new_env += v.b * w.b
            return new_env
        
        if isinstance(other, LinearEnvelope):
            return LinearEnvelope(lower = (self * other.lower).lower, upper=(self * other.upper).upper)


    #Suma de una expresión lineal con un escalar
    def __radd__(self, other): return self + other

    #Resta
    def __sub__(self, other): return self + (-other)

    #Resta Conmutativa
    def __rsub__(self, other): return -self + other

    #Negación
    def __neg__(self): return LinearExpression([-ai for ai in self.a], -self.b, self.x)

    #Multiplicación conmutativa
    def __rmul__(self, other): return self * other

    #División por escalar solamente
    def __truediv__(self, other): return self * other**(-1)

    #Transformación a cadena
    def __str__(self, constant=True) -> str:
        s = ""
        for i, ai in enumerate(self.a):
            if ai == 0: continue
            # Convert the representation from y to x
            if s!="": s+= " + "
            if ai == 1.0: s += f"x_{i+1}"
            else: s += f"{ai}*x_{i+1}"

        if constant:
            if s != "":
                s += f" + {str(self.b)}"
            else:
                s += str(self.b)
        return s

#Clase para Función Afín
class LinearEnvelope:
    def __init__(self, lower=None, upper=None, is_var=False, linear_expr=None):
        if linear_expr is not None:
            self.lower = linear_expr
            self.upper = linear_expr
        else:
            self.upper = upper
            self.lower = lower

        self.is_var = is_var

    def __rmul__(self, other):
        return None

    #Multiplicación
    def __mul__(self, other):
        if self.lower == None: return LinearEnvelope()
        if isinstance(other, (int, float, torch.Tensor)):
            if other >= 0:
                return LinearEnvelope(lower = self.lower * other, upper=self.upper * other)
            else:
                return LinearEnvelope(lower = self.upper * other, upper=self.lower * other)

        ## bilinear multiplication
        elif self.is_var and other.is_var:
            x1 = self.get_range()
            x2 = other.get_range()
            m2, m1 = x1.mid(), x2.mid() #derivadas parciales de x1*x2 en x_mid
            g = lambda x1, x2: x1 * x2 - m1 * x1 - m2 * x2
            g1 = g(x1.lb(), x2.lb())
            g2 = g(x1.lb(), x2.ub())
            g3 = g(x1.ub(), x2.lb())
            g4 = g(x1.ub(), x2.ub())
            g5 = None
            if m2 in x1 and m1 in x2:
                g5 = g(m2, m1) #derivada de g es 0 en este punto

            e1 = min((x for x in [g1, g2, g3, g4, g5] if x is not None), default=None)
            e2 = max((x for x in [g1, g2, g3, g4, g5] if x is not None), default=None)

            data_m1 = circular_queue.get_next()
            if data_m1:
                alfa_m1 = data_m1['alfa']
                data_m2 = circular_queue.get_next()
                alfa_m2 = data_m2['alfa']
            else:
                alfa_m1 = torch.tensor(1, dtype=torch.float64, requires_grad=True)
                alfa_m2 = torch.tensor(1, dtype=torch.float64, requires_grad=True)
                m_data = {
                    'alfa': alfa_m1,
                }
                m_data2 = {
                    'alfa': alfa_m2
                }
                circular_queue.add(m_data)
                circular_queue.add(m_data2)

            result = LinearEnvelope(lower = (m1*self + m2*other + e1).lower, upper = (m1*self + m2*other + e2).upper)
            return result

        # Multiplicacion de linear envelopes
        else:
            # no implementado retornar error
            x1 = self.get_range()
            if isinstance(other, LinearEnvelope):
                x2 = other.get_range()
                if x1.lb() >= 0 and x2.lb()>=0:
                    lower = (self.lower * other.lower).lower
                    upper = (self.upper * other.upper).upper
                    return LinearEnvelope(lower = lower, upper= upper)
            elif isinstance(other, LinearExpression):
                x2 = other()
                if x1.lb() >= 0 and x2.lb()>=0:
                    lower = (self.lower * other).lower
                    upper = (self.upper * other).upper
                    return LinearEnvelope(lower = lower, upper= upper)

            if (x1.lb() >= 0 and x2.ub() <= 0) or (x1.ub() <= 0 and x2.lb() >= 0) or (x1.ub() <= 0 and x2.ub() <= 0):
                v = self; w = other; sigpos=True
                if x1.ub() <= 0:
                    v = -self
                    sigpos = -sigpos

                if x2.ub() <= 0:
                    w = -other
                    sigpos = False

                if sigpos: 
                    return v*w
                else: 
                    return -(v*w)

            else: #evelopes que contienen el cero
                print(x1, x2)
                raise NotImplementedError

    #Negación que invierte el límite inferior y superior
    def __neg__(self):
        if self.lower == None: return LinearEnvelope()
        return LinearEnvelope(lower = -self.upper, upper=-self.lower)

    #Resta
    def __sub__(self, other):
        if self.lower == None: return -other
        return self + (-other)

    #Suma
    def __add__(self, other):
        if self.lower == None: return other
        #Con número
        elif isinstance(other, (int, float, torch.Tensor)):
            return LinearEnvelope(lower = self.lower + other, upper=self.upper + other)

        #Con intervalo
        elif isinstance(other, Interval):
            return LinearEnvelope(lower = self.lower + other.lb(), upper=self.upper + other.ub())

        #Con otra función afín
        elif isinstance(other, LinearEnvelope):
            return LinearEnvelope(lower = self.lower + other.lower, upper=self.upper + other.upper)

        elif isinstance(other, LinearExpression):
            return LinearEnvelope(lower = self.lower + other, upper=self.upper + other)

    #Exponente
    def __pow__(self, exponent):
        if self.lower == None: return LinearEnvelope()

        range = self.lower()
        lb, ub = range
        #Con un número (Solo 2 y 1)
        if isinstance(exponent, int):
            if exponent==2:
                return self.cheby_convex(f=lambda x: x**2, df_proj=lambda m: m/2, convex=True)
            elif exponent==-1:
                #para la división debemos revisar en qué cuadrante nos ubicamos para saber si la función es concava/convexa
                if lb > 0:
                    return self.cheby_convex(f=lambda x: 1/x, df_proj=lambda m: math.sqrt(-1/m), convex=True)
                elif ub < 0:
                    return self.cheby_convex(f=lambda x: 1/x, df_proj=lambda m: math.sqrt(-1/m), concave=True)
                else:
                    raise ValueError('LinearEnvelope**(-1) con intervalo que pasa por el 0.', self, range)
            elif lb > 0 and exponent == 3:
                return self.cheby_convex(f=lambda x: x**3, df_proj=lambda m: (m / 3) ** (1/2), convex=True)
            elif exponent == 4:
                return self.cheby_convex(f=lambda x: x**4, df_proj=lambda m: (m / 4) ** (1/3), convex=True)
            else:
                raise ValueError(f"exponent={exponent} is not implemented")


    def exp(self):
        return self.cheby_convex(
            f=lambda x: math.exp(x),
            df_proj=lambda m: math.log(m),
            convex=True)

    def sqrt(self):
        return self.cheby_convex(
            f=lambda x: math.sqrt(x),
            df_proj=lambda m: 1 / (4 * math.pow(m, 2)),
            concave=True)


    #División
    def __truediv__(self, other):
        if self.lower == None: return LinearEnvelope()
        return self * other**(-1)

    def __rtruediv__(self, other):
        if self.lower == None: return LinearEnvelope()
        return other * self**(-1)

    #Aproximación de Chebyshev para funciones convexas y cóncavas
    def cheby_convex(self, f=lambda x: x**2, df_proj=lambda m: m/2, convex=False, concave=False):
        x = self.get_range()
        lb, ub = x.lb(), x.ub()

        # Determinar la Pendiente m
        m = (f(ub) - f(lb)) / (ub - lb)
        m_static = m

        # Identificar el Punto de Tangencia xp
        # Necesitamos resolver df(xp) = m
        xp = df_proj(m)
        if xp > ub: xp = ub
        elif xp < lb: xp = lb

        if convex: 
            if m >= (f(ub) - f(lb)) / (ub - lb):
                upper = m_static*self.upper + (f(lb) - m_static * lb)
            else:
                upper = m_static*self.upper + (f(ub) - m_static * ub)
            lower = m*self.lower + (f(xp) - m * xp)
        elif concave:
            if m >= (f(ub) - f(lb)) / (ub - lb):
                lower = m_static*self.lower + (f(ub) - m_static * ub)
            else:
                lower = m_static*self.lower + (f(lb) - m_static * lb)

            upper = m*self.upper + (f(xp) - m * xp)

        return LinearEnvelope(lower = lower, upper= upper)

    #Multiplicación conmutativa
    def __rmul__(self, other):
        if self.lower == None: return LinearEnvelope()
        return self * other

    #Conversión a intervalo
    def get_range(self):
        return IntervalTensor(self.lower(Affine.x_ref).lb(), self.upper(Affine.x_ref).ub())

    #Pasar a cadena para mostrar
    def __str__(self):
        if not self.upper:
            return "Not self.lower"
        if not self.lower or not self.upper:
            return ""
        if self.lower.a == self.upper.a:
            return f"{self.lower.__str__(constant=False)} + [{self.lower.b}, {self.upper.b}]"
        else:
            return f"[{self.lower}, {self.upper}]"