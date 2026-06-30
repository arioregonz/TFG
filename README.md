# TFG

# Plataforma de simulación para la caracterización acústica de la madera mediante ultrasonidos

Trabajo Fin de Grado — Grado en Ingeniería de Sistemas Audiovisuales y Multimedia
Universidad Rey Juan Carlos (URJC)

## Descripción

Plataforma de simulación de la propagación de ondas ultrasónicas en madera, basada
en el Método de los Elementos Espectrales (SEM). Permite modelar secciones
transversales de troncos de distintas especies, con y sin defectos internos, y
evaluar la detección y localización de dichos defectos mediante el tiempo de vuelo
(TOF) de la señal.

El trabajo reproduce y amplía los resultados de Espinosa et al. (2019), añadiendo
el estudio del tipo de relleno del defecto (aire frente a agua), la comparación
entre especies (roble rojo y pino de hoja larga) y la localización angular del
defecto mediante múltiples receptores.

## Requisitos

- **SPECFEM2D** — solver de elementos espectrales (propagación de ondas)
- **Gmsh** (≥ 4.4) — generación de mallas externas
- **MATLAB** — postprocesado de sismogramas y generación de figuras
- **Python 3** (con NumPy) — generación de la señal de excitación

## Estructura del repositorio

```
.
├── senal/              Generación de la señal chirp (Python)
├── mallas/             Geometrías de Gmsh (.geo)
├── configuracion/      Ficheros Par_file y señal de excitación
├── scripts_matlab/     Postprocesado y validación (MATLAB)
└── stations/           Configuración de receptores
```

## Cómo reproducir los resultados

1. Generar la señal de excitación:
   ```
   python3 senal/generar_chirp.py
   ```

2. Generar la malla con Gmsh (geometría circular):
   ```
   gmsh mallas/tronco_circular.geo -2 -format msh2 -o tronco.msh
   ```

3. Ejecutar la simulación con SPECFEM2D (ver Par_file en configuracion/).

4. Procesar los resultados en MATLAB:
   ```
   run('scripts_matlab/validacion_base.m')
   ```

## Parámetros de las especies

| Parámetro | Roble rojo | Pino de hoja larga |
|-----------|------------|--------------------|
| ρ (kg/m³) | 706 | 661 |
| E_R (GPa) | 2.54 | 1.537 |
| E_T (GPa) | 1.35 | 0.829 |
| V radial (m/s) | 1897 | 1525 |
| V tangencial (m/s) | 1383 | 1120 |

## Autoría

Autora: Ariana Marjorie Ore Gonzales
Tutor: Víctor Francisco Martín Martínez
Curso 2025-2026

## Licencia

MIT License
