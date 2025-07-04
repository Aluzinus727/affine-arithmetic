from solvers.af1 import AF1
from solvers.af2 import AF2
from solvers.variant import LinearExpression
from solvers.ibex import get_optimal
from utils import process_data, export_to_excel, str_to_func
from instances import instances

evaluators = [
    ("RevisedAF", LinearExpression),
    ("AF1", AF1),
    ("AF2", AF2),
]

results = []
new_results = {}

for instance in instances:
    fname, f, intervals = instance
    new_results[fname] = {}

    for interval in intervals:
        print('Testing con', fname, interval)
        interval_str = ''.join((str(x) for x in interval))
        new_results[fname][interval_str] = {}

        # ---------- IBEX ----------
        optimal_linear, optimal_range = get_optimal(f, interval)
        ibex_data = {
            "Function": f,
            "Interval": interval_str,
            "Implementation": "Ibex",
            "Min": min(optimal_range),
            "Max": max(optimal_range)
        }
        new_results[fname][interval_str]['Ibex'] = ibex_data
        results.append(ibex_data)

        for ev_name, EvClass in evaluators:
            print('Ev_name', ev_name) 
            converted_vars = EvClass.to_intervals(interval) # Intervalos como [1, 2] a LinearExpression, AF1, AF2...
            print('converted vars', converted_vars)
            f_modified = f.replace('[', '').replace(']', '').replace('^', '**') # Func como string. Formato sympy
            print('f_modified', f_modified)
            func = str_to_func(f_modified, len(converted_vars)) # Callable
            result = func(*converted_vars)
            f_range = result.get_range()

            print('result, f_range, ev_name', result, f_range, ev_name)

            ev_data = {
                "Function": f,
                "Interval": interval_str,
                "Implementation": ev_name,
                "Min": min(f_range),
                "Max": max(f_range)
            }

            new_results[fname][interval_str][ev_name] = ev_data
            results.append(ev_data)

            EvClass.reset()

raw_data, _, error_summary, method_stats = process_data(new_results)
export_to_excel(raw_data, 'raw.xlsx')
export_to_excel(error_summary, 'error_summary.xlsx')
export_to_excel(method_stats, 'method_stats.xlsx')
# detailed_rows, function_summary_rows, method_overall_rows, raw_rows = process_data(new_results)
# export_to_excel(detailed_rows, 'details.xlsx')
# export_to_excel(function_summary_rows, 'summary_by_function.xlsx')
# export_to_excel(method_overall_rows, 'summary_by_method.xlsx')
# export_to_excel(raw_rows, 'raw.xlsx')