from pyibex import Interval, IntervalVector
from typing import List

class AFPaper:
    def __init__(
        self, 
        a: List, 
        x: List, 
        error: Interval = Interval(0)
    ):
        self.a = a
        self.x = x

        self.error = error

    def to_interval(self):
        sum_terms = 0
        for idx, a in enumerate(self.a):
            sum_terms += a * self.x[idx]

        return sum_terms + self.error

    def __add__(self, other):
        if isinstance(other, AFPaper):
            new_a = [a + other.a[idx] for idx, a in enumerate(self.a)]
            new_x = self.x
            new_error = self.error + other.error

            return AFPaper(
                a=new_a, 
                x=new_x,
                error=new_error
            )

        elif isinstance(other, (int, float)):
            new_a = self.a
            new_x = self.x
            new_error = self.error + other

            return AFPaper(
                a=new_a, 
                x=new_x,
                error=new_error
            )
        elif isinstance(other, Interval):
            new_a = self.a
            new_x = self.x
            new_error = self.error + other

            return AFPaper(
                a=new_a, 
                x=new_x,
                error=new_error
            )

    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        if isinstance(other, AFPaper):
            new_a = [a - b for (a, _), (b, _) in zip(self.a, other.a)]
            new_x = self.x
            new_error = self.error - other.error

            return AFPaper(
                a=new_a, 
                x=new_x,
                error=new_error
            )

    def __rsub__(self, other):
        if isinstance(other, AFPaper):
            return other.__sub__(self)

    def __mul__(self, other):
        if isinstance(other, (int, float, Interval)):
            new_a = [other * a for a in self.a]
            new_x = self.x
            new_error = other * self.error

            return AFPaper(
                a=new_a, 
                x=new_x,
                error=new_error
            )

        elif isinstance(other, AFPaper):
            sum_1 = 0
            for idx, a in enumerate(self.a):
                quad_term_a = [1 if i == idx else 0 for i in range(len(self.a))]
                quad_term_af = AFPaper(a = quad_term_a, x = self.x, error = Interval(0))
                quad_term = quad_term_af**2 
                sum_1 += a * other.a[idx] * quad_term

            sum_2 = 0
            for idx, a in enumerate(self.a):
                for jdx, b in enumerate(other.a):
                    if idx == jdx:
                        continue

                    x_i = self.x[idx]
                    x_j = self.x[jdx]
                    bilinear_term = self.bilinear_multiplication(x_i, x_j)
                    sum_2 += a * b * bilinear_term

            # Error e_u
            m_v = self.error.mid()
            m_w = other.error.mid()

            e_v = m_v + (self.error - m_v)
            e_w = m_w + (other.error - m_w)

            # First two terms
            sum_terms = 0
            terms = []
            for idx, a in enumerate(self.a):
                b = other.a[idx]

                # Según la fórmula, en cada iter se debería multiplicar por x, pero el resultado daría un intervalo
                # Se calculan los terminos y se transforman a AF para operar.
                # De esta forma, queda la multiplicación implicita dentro de la instancia:
                # a_1 = (m_w * a_1) + (m_v * b_1) -> a_1 * x_1
                new_term = (m_w * a) + (m_v * b)
                terms.append(new_term)
            
            sum_terms = AFPaper(a=terms, x=self.x, error=Interval(0))

            e_u = (e_v - m_v) * other.to_interval() + (e_w - m_w) * self.to_interval() + m_v * e_w + m_w * e_v

            # Para debug y ver cada uno de los 4 términos generados
            # print('Sum de i=1 hasta n (m_w * a_i + m_v * b_i) * x_i', sum_terms)
            # print('Sum de i=1 hasta n (a_i*b_i*(x_i)^2)', sum_1)
            # print('Doble sum', sum_2)
            # print('eu: ', e_u)

            result = sum_terms + sum_1 + sum_2 + e_u
            return result

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return self.__mul__(other)

    def bilinear_multiplication(self, x1: Interval, x2: Interval):
        if isinstance(x1, Interval) and isinstance(x2, Interval):
            m1 = x2.mid()
            m2 = x1.mid()

            g = lambda x_1, x_2: x_1 * x_2  - (m1 * x_1 + m2 * x_2)

            # print('m1', m1)
            # print('m2', m2)
            # print(g(x1.lb(), x2.lb()), g(x1.lb(), x2.ub()), g(x1.ub(), x2.lb()), g(x1.ub(), x2.ub()), -m1*m2)

            eu_lb = min(
                g(x1.lb(), x2.lb()),
                g(x1.lb(), x2.ub()),
                g(x1.ub(), x2.lb()),
                g(x1.ub(), x2.ub()),
                -m1*m2
            ) 

            # print('eu_lb', eu_lb)

            eu_ub = max(
                g(x1.lb(), x2.lb()),
                g(x1.lb(), x2.ub()),
                g(x1.ub(), x2.lb()),
                g(x1.ub(), x2.ub()),
                -m1*m2
            ) 

            new_a = [m1, m2]
            new_x = [x1, x2]
            new_error = Interval(eu_lb, eu_ub)

            return AFPaper(
                a = new_a,
                x = new_x,
                error = new_error
            )

    def non_affine(self, f, df_proj):
        interval = self.to_interval()
        lb = interval.lb()
        ub = interval.ub()

        m = (f(ub) - f(lb)) / (ub - lb)
        xp = df_proj(m)

        g = lambda v: f(v) - m * v

        eu_lb = min(g(xp), g(lb), g(ub))
        eu_ub = max(g(xp), g(lb), g(ub))

        return m * self + Interval(eu_lb, eu_ub)

    def __pow__(self, exponent):
        if exponent == 2:
            return self.non_affine(f=lambda x: x**2, df_proj=lambda m: m/2)
        elif exponent == 3:
            return self.cheby_convex(f=lambda x: x**3, df_proj=lambda m: math.sqrt(m/3))
        elif exponent == 4:
            return self.cheby_convex(f=lambda x: x**4, df_proj=lambda m: (m / 4) ** (1 / 3))
        elif exponent == 6:
            return self.cheby_convex(f=lambda x: x**6, df_proj=lambda m: (m / 6) ** (1 / 5))
        elif exponent == 7:
            return self.cheby_convex(f=lambda x: x**7, df_proj=lambda m: (m / 7) ** (1 / 6))


    def __str__(self):
        lb, ub = self.error.lb(), self.error.ub()
        return " + ".join(f"{coef} * x{i+1}" for i, coef in enumerate(self.a)) + f" + [{lb}, {ub}]"



