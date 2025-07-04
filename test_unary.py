from solvers.variant import LinearExpression
from gradient_descent import grad_descent
from instances import instances
from utils import str_to_func

instances_unary = [
    [    
        "raiz cuadrada",
        "sqrt(x0)",
        [
           [[1, 2]]
        ]
    ],
    [
        "cuadrado",
        "x0**2",
        [
            [[1, 2]]
        ],
    ],
    [
        "exponencial",
        "exp(x0)",
        [
            [[1, 2]]
        ]
    ],
    [
        "cos",
        "cos(x0)",
        [
            [[2, 2.25]]
        ]
    ],
]


for fname, fstr, test_cases in instances[1:]:
    print(f'Probando funcion {fname}, {fstr}')
    for test in test_cases:
        fstr.replace('[', '').replace(']', '').replace('^', '**')
        func = str_to_func(fstr, len(test))

        result = grad_descent(test, func)
        f_range = result.get_range()

        print(f_range)

    break

    

for fname, fstr, test_cases in instances:
    for test in test_cases:
        fstr.replace('[', '').replace(']', '').replace('^', '**')
        func = str_to_func(fstr, len(test))

        converted_vars = LinearExpression.to_intervals(test)

        result = func(*converted_vars)
        # print(type(result))
        # print('result nuevo', result)
        f_range = result.get_range()

        print(f_range)




