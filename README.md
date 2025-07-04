# AF Revisada
## Instalación
IBEX no corre en Windows, por lo que para correr esto en Windows hay que utilizar WSL. Tutorial cortesía del Eduardo
1. Instalar WSL (esto automáticamente instala Ubuntu):  
- Abrir Símbolo del sistema como administrador, y ejecutar el siguiente comando:  
_wsl --install_
2. Abrir WSL, y ejecutar los siguientes comandos:  
_sudo apt update && sudo apt upgrade_
Acá el tutorial dice que hay que instalar python con sudo, pero esto a mi **no** me gusto funcionó. Hay que fijarse con que versión de python fue compilado ibexopt. En mi caso, lo que funcionó fue instalar la versión 3.11.6 con [PyEnv](https://github.com/pyenv/pyenv).
3. Abrir WSL y ejecutar los siguientes comandos:  
_wget -O pybex1.9.0.tar.gz (URL de google drive con la versión deseada)
pip install pybex pybex1.9.0.tar.gz &> /dev/null_

## Solvers:

Del paper Extensions of Affine Arithmetic: Application to Unconstrained Global Optimization:
- solvers/af1.py implementa AF1
- solvers/af2.py implementa AF2

De otros papers:
- solvers/variant implementa AF Revisada con LinearExpression y LinearEnvelope, según el proyecto de título
- solvers/paper.py implementa (a medias) AF revisada según el paper que se está escribiendo. Falta agregar las funciones sqrt, exp, log, sin, cos, pero es cosa de verificar que el intervalo sea concavo o convexo en cada funcióny pasarle f(x) y df_proj a la función non_affine.

## Instancias
Todas las instancias (utilizadas y no utilizadas) están en instances.py. Cada instancia tiene 5 casos de prueba con intervalos elegidos específicamente para que las operaciones no afines sean siempre sobre variables estrictamente cóncavas o convexas.

 ## Experimentos
- gradient_decent.py ocupa la implementación de AF Revisada con LinearExpression y LinearEnvelope para generar las expresiones iniciales y aplica gradiente descendente. Importa la variante desde solvers/variant_gd.py, que en principio es lo mismo que está en solvers/variant.py pero ocupa tensores para luego optimizarlos con el gradiente descendente.
- paper_af_rev.py ocupa la implementación de AF Revisada del paper que se está escribiendo y la compara con la implementación de AF2. Es el ejemplo de multiplicación que se puso en el paper.
- paper_af2.py ocupa un ejemplo del paper desde donde se obtuvo la AF1 y AF2 para comparar. Este ejemplo no da el mismo resultado del paper, falta revisar.
- test_unary.py prueba las funciones no afines básicas en ambas implementaciones de AF revisada (la que funciona con tensores y la normal)
- exp.py tiene los experimentos completos (instancias de 2 a 10 variables). Ocupa AF1, AF2, Ibex y AFRevised para comparar y los exporta a excel, agrupados por función y método.

## Problemas
1. El gradiente descendente funciona, pero solo probé con las funciones sin multiplicaciones no afines que los rangos resultantes fueran válidos, no probé con los que tiene multiplicación. 
Si el código del gradiente descendente da problemas y tira error de variables fuera de dominio o cosas parecidas, probar bajar el learning rate antes de debugear, algunas funciones necesitan un lr bajísimo para funcionar y no salirse de dominio (5e-10 o cercano). Evidentemente, la mejora es mínima, pero al menos se nota que no es el código el malo y se nota que va mejorando el resultado.
2. La AF2 no da los resultados esperados (en paper_af2.py), se supone que al agregar en AF2 se agregan los términos de error negativos [-1, 0] y [0, 1] para intentar solucionar problemas de sobreestimación en potencias pares/positivas, ya que los métodos tradicionales siempre incluyen números negativos en sus soluciones, por lo que el intervalo se puede acotar en teoría porque las soluciones siempre serán positivas. Sin embargo, en el paper no especifica como lo hace, y de las cosas que probé, si bien algunas me daban resultados similares, la expresión que se genera no es la misma, hay que revisar la implementación de \_\_pow\_\_.
3. Hay que revisar la función de exportar a excel para asegurarse que calcula bien los porcentajes de error y la desviación estándar. O sino, exportar el raw data y hacer los cálculos manuales en excel.
