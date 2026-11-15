# Proyecto de Sistemas de Control

Proyecto académico de identificación y control de velocidad de un motor de corriente directa, desarrollado en MATLAB y Simulink para el curso de Sistemas de Control de la Universidad de Costa Rica.

Incluye datos experimentales, identificación de la planta, comparación de controladores PI por síntesis analítica, Chien y Fertik-Sharpe, y un controlador mediante lugar geométrico de las raíces (LGR) discreto. Los análisis consideran seguimiento de referencia, rechazo de perturbaciones, desempeño y robustez.

## Requisitos

- MATLAB y Simulink. El modelo incluido fue guardado con R2025b; utiliza esa versión o una posterior para abrirlo.
- Control System Toolbox para funciones de transferencia y análisis de sistemas.
- System Identification Toolbox para cargar y trabajar con el modelo identificado y la sesión de identificación.

## Archivos principales

Los scripts, el modelo y los datos se encuentran juntos en `Proyecto Control Grupo 3-4 solo matlab/Proyecto Control Grupo 3-4/`.

| Archivo | Contenido |
| --- | --- |
| `PruebaPlantas.m` | Identificación de la planta a partir de `DatosMotor.csv`, comparación de modelos y métricas de ajuste. |
| `P1D.mat` | Modelo identificado utilizado por `PruebaPlantas.m`. |
| `mydata.sid` | Sesión de System Identification. |
| `proyectoControl.slx` | Simulación de los controladores continuos y del controlador discreto. |
| `proyectoControlcontrolador.m` | Gráficas y métricas de las señales obtenidas de la simulación. |
| `curvas_expn.m` | Comparación de las respuestas experimentales de los tres controladores PI. |
| `ms.m` | Análisis de sensibilidad máxima, márgenes y diagramas de Bode. |
| `bodeProyecto.m` | Diagramas de Bode de la inversa del lazo abierto para los controladores PI. |
| `*.csv` | Mediciones de identificación y de las pruebas experimentales. |

El [informe del proyecto](Proyecto%20Control%20Grupo%203-4%20solo%20matlab/ProyectoControl3-4.pdf) contiene el desarrollo y los resultados. También se incluye la imagen de referencia `fnctransferenciatoolbox.png` en la raíz del repositorio.

## Uso

Clona el repositorio o descarga su contenido. En MATLAB, sitúate en la raíz del repositorio y cambia a la carpeta de trabajo:

```matlab
cd(fullfile('Proyecto Control Grupo 3-4 solo matlab', 'Proyecto Control Grupo 3-4'))
```

Para identificar y comparar los modelos de la planta:

```matlab
PruebaPlantas
```

Para simular los controladores y analizar las señales resultantes:

```matlab
out = sim('proyectoControl');
proyectoControlcontrolador
```

El script `proyectoControlcontrolador.m` necesita la variable `out` generada por la simulación. Para inspeccionar el diagrama, ejecuta `open_system('proyectoControl')`.

Los análisis experimentales y de frecuencia pueden ejecutarse por separado:

```matlab
curvas_expn
ms
bodeProyecto
```

Ejecuta cada análisis desde la carpeta de trabajo para que encuentre los archivos de datos. Algunos scripts limpian el espacio de trabajo y cierran figuras; ejecuta cada uno por separado si quieres conservar sus resultados.

## Autores

- Walther Barrantes Yan
- Kenneth Vizcaíno Jiménez
- Andrés Víquez Zamora

Escuela de Ingeniería Eléctrica, Universidad de Costa Rica.

## Licencia

Este proyecto se distribuye bajo la [licencia MIT](LICENSE).
