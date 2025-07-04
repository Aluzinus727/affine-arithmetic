from pyibex import Interval, IntervalVector
from solvers.variant import LinearExpression
from utils import str_to_func, init_boxes
from memory import Affine
import os, re, subprocess

def bounds(vars, xl, x, ctrs=None, prec=1e-4):
    #minimization
    # Abre el archivo en modo escritura
    with open('tmp.txt', 'w') as f:
        f.write('variables\n')
        f.write(vars+' in '+str(x)+';\n\n')
        f.write('minimize '+xl+';')
        if ctrs is not None:
            f.write('\nconstraints\n')
            for ctr in ctrs: f.write(ctr+'\n')
            f.write("end\n")


    p = subprocess.Popen(["./ibexopt", '-a', str(prec), '-r', str(prec), '-t', '5', 'tmp.txt'], stdout=subprocess.PIPE)

    for line in p.communicate()[0].decode().splitlines():
        if line.find("f* in")!=-1:
            lb = float(re.search('\[(.*),(.*)\]', line).groups(0)[1])
        if line.find("number of cells:")!=-1:
            n_cells = int(re.search('\t\t([0-9]*)', line).groups(0)[0])

    #maximization
    with open('tmp.txt', 'w') as f:
        f.write('variables\n')
        f.write(vars+' in '+str(x)+';\n\n')
        f.write('minimize -('+xl+');\n')
        if ctrs is not None:
            f.write('\nconstraints\n')
            for ctr in ctrs: f.write(ctr+'\n')
            f.write("end\n")


    p = subprocess.Popen(["./ibexopt", '-a', str(prec), '-r', str(prec), '-t', '5', 'tmp.txt'], stdout=subprocess.PIPE)

    for line in p.communicate()[0].decode().splitlines():
        if line.find("f* in")!=-1:
            ub = -float(re.search('\[(.*),(.*)\]', line).groups(0)[1])
        if line.find("number of cells:")!=-1:
            n_cells += int(re.search('\t\t([0-9]*)', line).groups(0)[0])


    return (lb,ub), n_cells


def optimal_linearization(x_ref=IntervalVector([Interval(0,2),Interval(0,2)]), real_function="x_1^2+x_2^2", linear_aprox="2*x_1+2*x_2", error=None):
    x=f"{x_ref}" #.replace("(","[").replace(";",",").replace(")","]")

    real_function = re.sub(r'x_(\d+)', lambda m: f'x[{int(m.group(1)) - 1}]', real_function)
    linear_aprox_ = re.sub(r'x_(\d+)', lambda m: f'x[{int(m.group(1)) - 1}]', linear_aprox)

    intercepts,_ = bounds(f"x[{x_ref.size()}]", f"{real_function}-({linear_aprox_})", x, prec=1e-8)

    if error is None:
        return f"{linear_aprox}+[{intercepts[0]},{intercepts[1]}]"
    elif error == "lb":
        return f"{linear_aprox}+{intercepts[0]}"
    elif error == "ub":
        return f"{linear_aprox}+{intercepts[1]}"


def get_lower_bound(expr: str, x_values: dict[str, tuple[float, float]]) -> float:
    # 1. Extraer el intervalo del final
    interval_match = re.search(r"\[ *(-?\d+(?:\.\d+)?), *(-?\d+(?:\.\d+)?) *\]\s*$", expr)
    if not interval_match:
        print("⚠ No se detectó intervalo. Expr fue:", repr(expr))
        raise ValueError("No se encontró un intervalo al final de la expresión.")
    interval_lb = float(interval_match.group(1))

    # 2. Quitar el intervalo del string para quedarnos solo con la expresión simbólica
    expr_core = expr[:interval_match.start()].strip()

    # 3. Buscar todos los coeficientes y variables, como 3.0*x_1 o -2*x_3
    terms = re.findall(r"([+-]?\s*\d*\.?\d*)\s*\*\s*(x_\d+)", expr_core)
    
    total = 0.0
    for coef_str, var in terms:
        coef = float(coef_str.replace(" ", "") or "1")
        var_min = x_values[var][0]
        total += coef * var_min

    # 4. Agregar el extremo inferior del intervalo
    return total + interval_lb

def get_upper_bound(expr: str, x_values: dict[str, tuple[float, float]]) -> float:
    interval_match = re.search(r"\[ *(-?\d+(?:\.\d+)?), *(-?\d+(?:\.\d+)?) *\]\s*$", expr)
    if not interval_match:
        print("⚠ No se detectó intervalo. Expr fue:", repr(expr))
        raise ValueError("No se encontró un intervalo al final de la expresión.")
    interval_ub = float(interval_match.group(2))

    expr_core = expr[:interval_match.start()].strip()
    terms = re.findall(r"([+-]?\s*\d*\.?\d*)\s*\*\s*(x_\d+)", expr_core)

    total = 0.0
    for coef_str, var in terms:
        coef = float(coef_str.replace(" ", "") or "1")
        var_max = x_values[var][1]
        total += coef * var_max

    return total + interval_ub



def get_optimal(f, intervals):
    glob = []

    vars = [Interval(min(inter[0], inter[1]), max(inter[0], inter[1])) for inter in intervals]
    new_x = IntervalVector(vars)

    x, var_envs = init_boxes(new_x, vars)
    func = str_to_func(f, len(intervals))
    experimental_res = func(*x)
    optimal_res_lower = optimal_linearization(new_x, f.replace('**', '^'), experimental_res.lower.__str__(constant=False))

    eval_vars = {}
    for i, interval in enumerate(intervals, 1):
        eval_vars[f"x_{i}"] = interval

    lb = get_lower_bound(optimal_res_lower, eval_vars)
    ub = get_upper_bound(optimal_res_lower, eval_vars)
    
    return optimal_res_lower, (lb, ub)


