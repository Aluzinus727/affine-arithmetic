import torch
from solvers.variant_gd import LinearExpression
from instances import instances
from utils import str_to_func, circular_queue

def grad_descent(boxes, func, epochs=25000, lr=5e-10):
    """
    Función principal del gradiente descendiente
    Cada función funciona bien con un lr distinto, ir probando
    Hay algunas que funcionan bien con lrs más altos (5e-2, 5e-5)
    Otras funcionan bien solo con lrs bajísimos (5e-9, 5e-11)
    """
    torch.autograd.set_detect_anomaly(True)
    # print(f'Evaluando con intervalos {boxes}, func {func}')
    intervals = LinearExpression.to_intervals(boxes)

    circular_queue.reset()
    z_af = func(*intervals)
    circular_queue.full = True
    all = circular_queue.get_all()
    z_int = z_af.get_range()
    print(f'Expresión afín inicial: {z_af}. LB = {z_int.lower} UB = {z_int.upper}')

    params = []
    items = circular_queue.get_all()
    for item in items:
        alfa = item['alfa']
        if alfa:
            params.append(alfa)

    optimizer = torch.optim.SGD(params, lr=lr, maximize=False)
    old_loss = float('inf')

    # Loop principal del gradiente descendente
    for epoch in range(epochs):
        loss = z_int.upper - z_int.lower

        if epoch % 100 == 0:
            print(f'[EPOCH {epoch}]. LB = {z_int.lower} UB = {z_int.upper}. Loss = {loss}')

        if loss > old_loss:
            # print("Loss increased. Stopping early at epoch", epoch)
            break  # Termina la iteración si la pérdida actual es mayor que la anterior

        optimizer.zero_grad()
        loss.backward(retain_graph=True)
        optimizer.step() # Paso

        old_loss = loss

        z_af = func(*intervals)
        z_int = z_af.get_range()

    return z_af

for fname, fstr, test_cases in instances:
    print(f'Probando funcion {fname}, {fstr}')
    for test in test_cases:
        fstr.replace('[', '').replace(']', '').replace('^', '**')
        func = str_to_func(fstr, len(test))

        # try:
        result = grad_descent(test, func)
        f_range = result.get_range()

        print(f_range)

