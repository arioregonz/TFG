# Plataforma de simulación para la caracterización acústica de la madera mediante ultrasonidos

Trabajo Fin de Grado — Grado en Ingeniería de Sistemas Audiovisuales y Multimedia
Universidad Rey Juan Carlos (URJC), ETSI de Fuenlabrada

Plataforma de simulación de la propagación de ondas ultrasónicas en secciones
transversales de troncos, basada en el Método de los Elementos Espectrales (SEM)
mediante SPECFEM2D y Gmsh. Reproduce y amplía los resultados de Espinosa et al.
(2019), incorporando el estudio del tipo de relleno del defecto (aire frente a
agua), la comparación entre especies (roble rojo y pino de hoja larga) y la
localización angular del defecto con múltiples receptores.

---

## Requisitos

- **SPECFEM2D** — solver de elementos espectrales (propagación de ondas)
- **Gmsh** (≥ 4.4) — generación de mallas externas
- **MATLAB** — postprocesado de sismogramas y generación de figuras
- **Python 3** (con NumPy) — generación de la señal de excitación y utilidades

---

## Estructura general de carpetas

El proyecto se organiza en dos especies (`Roble/` y `Pino/`), cada una con dos
estrategias de mallado:

```
EXAMPLES/
├── Roble/
│   ├── Mesh_interna/        Mallas generadas por el mallador interno de SPECFEM2D
│   └── Mesh_externa/        Mallas generadas con Gmsh (geometría circular y rectangular)
└── Pino/
    ├── Mesh_interna/
    └── Mesh_externa/
```

- **`Mesh_interna/`** contiene los modelos de tronco **rectangular** mallados con
  el mallador interno de SPECFEM2D. Se emplean para la validación física y el
  estudio paramétrico de tamaños y posiciones del defecto.
- **`Mesh_externa/`** contiene los modelos mallados con **Gmsh**, incluyendo la
  geometría **circular** (tronco realista) y la geometría rectangular con PML
  transfinito para la validación cruzada de mallas.

---

## Carpetas por caso de estudio (Roble)

La misma estructura se replica para el pino sustituyendo `Roble_` por `Pino_`.

### Malla interna — validación y estudio paramétrico

| Carpeta | Caso de estudio |
|---|---|
| `Mesh_interna/Roble_sin_hueco_Norm` | Tronco sano (referencia de validación física) |
| `Mesh_interna/Roble_con_hueco_N_C` | Hueco centrado 10×10 cm (validación vs. tesis, Fig.63) |
| `Mesh_interna/Roble_con_hueco_N_C_Mediano` | Hueco centrado 5×5 cm (barrido de tamaños) |
| `Mesh_interna/Roble_con_hueco_N_C_Pequeño` | Hueco centrado 2×2 cm (barrido de tamaños) |
| `Mesh_interna/Roble_con_hueco_N_UR` | Hueco 10×10 cm arriba-derecha (barrido de posiciones) |
| `Mesh_interna/Roble_con_hueco_N_UL` | Hueco 10×10 cm arriba-izquierda (barrido de posiciones) |
| `Mesh_interna/Roble_con_hueco_N_DL` | Hueco 10×10 cm abajo-izquierda (barrido de posiciones) |
| `Mesh_interna/Roble_con_hueco_N_DR` | Hueco 10×10 cm abajo-derecha (barrido de posiciones) |

### Malla externa — geometría rectangular (validación de mallas)

| Carpeta | Caso de estudio |
|---|---|
| `Mesh_externa/Roble_R_sin_hueco_PML_T` | Tronco rectangular sano con PML transfinito (validación malla interna vs. externa) |

### Malla externa — geometría circular (tronco realista)

| Carpeta | Caso de estudio |
|---|---|
| `Mesh_externa/Roble_C_sin_hueco_N_2` | Tronco circular sano (validación geometría realista) |
| `Mesh_externa/Roble_C_con_hueco_N_3` | Hueco cilíndrico centrado Ø10 cm, relleno de **aire** |
| `Mesh_externa/Roble_C_hueco_agua` | Hueco centrado Ø10 cm, relleno de **agua** (pudrición húmeda) |
| `Mesh_externa/Roble_C_hueco_5cm` | Hueco cilíndrico centrado Ø5 cm (barrido de tamaños) |
| `Mesh_externa/Roble_C_hueco_2cm` | Hueco cilíndrico centrado Ø2 cm (barrido de tamaños) |
| `Mesh_externa/Roble_C_hueco_UR` | Hueco Ø10 cm arriba-derecha (barrido de posiciones) |
| `Mesh_externa/Roble_C_hueco_UL` | Hueco Ø10 cm arriba-izquierda (barrido de posiciones) |
| `Mesh_externa/Roble_C_hueco_DL` | Hueco Ø10 cm abajo-izquierda (barrido de posiciones) |
| `Mesh_externa/Roble_C_hueco_DR` | Hueco Ø10 cm abajo-derecha (barrido de posiciones) |

### Malla externa — localización con 16 receptores

| Carpeta | Caso de estudio |
|---|---|
| `Mesh_externa/Roble_C_16sensores_sano` | Corona de 16 receptores, tronco sano (referencia) |
| `Mesh_externa/Roble_C_16sensores_hueco` | 16 receptores, hueco centrado (patrón simétrico) |
| `Mesh_externa/Roble_C_16sensores_exc` | 16 receptores, hueco excéntrico (localización angular) |

### Estructura interna de cada carpeta de caso

```
Roble_C_con_hueco_N_3/
├── DATA/
│   ├── Par_file              Configuración de la simulación
│   ├── SOURCE                Definición de la fuente (posición y anglesource)
│   ├── STATIONS              Definición de receptores (solo casos multisensor)
│   ├── chirp_source.txt      Señal de excitación (chirp externo)
│   ├── three.geo             Geometría de Gmsh (solo malla externa)
│   └── three.msh             Malla generada (solo malla externa)
├── OUTPUT_FILES/
│   ├── AA.S0001.BXX.semd     Sismograma del receptor (componente X)
│   └── forward_image*.jpg    Instantáneas del frente de onda (snapshots)
└── run_this_example.sh       Script de ejecución de la simulación
```

---

## Parámetros de las especies

Constantes elásticas empleadas en el `Par_file` (material anisótropo tipo 2).

| Parámetro | Roble rojo | Pino de hoja larga |
|---|---|---|
| Densidad ρ (kg/m³) | 706 | 661 |
| c11 = E_R (GPa) | 2.54 | 1.537 |
| c33 = E_T (GPa) | 1.35 | 0.829 |
| c55 = G_RT (GPa) | 0.319 | 0.181 |
| c13 (GPa) | 0.758 | 0.429 |
| V radial (m/s) | 1897 | 1525 |
| V tangencial (m/s) | 1383 | 1120 |

Líneas del `Par_file` (formato: `número tipo rho c11 c13 c15 c33 c35 c55 c12 c23 c25 0 0 0`):

```
# Roble rojo:
2 2 706.d0 2.54d9 7.5832d8 0 1.35d9 0 3.19d8 7.58d8 7.58d8 0 0 0 0
# Pino de hoja larga:
2 2 661.d0 1.537d9 4.289d8 0 8.29d8 0 1.81d8 4.289d8 4.289d8 0 0 0 0
```

Materiales de relleno del defecto (material acústico tipo 1):

```
# Aire:  1 1 1.2d0   343.d0  0 0 0 9999 9999 0 0 0 0 0 0
# Agua:  1 1 1000.d0 1480.d0 0 0 0 9999 9999 0 0 0 0 0 0
```

---

## Comandos

Esta sección recoge, paso a paso, los comandos concretos para reproducir un caso
completo de la malla externa circular, explicando el porqué de cada uno.

### 1. Generar la señal de excitación (chirp)

La fuente emplea una señal chirp con ventana gaussiana (22–50 kHz, 45 µs),
generada en Python y leída por SPECFEM2D como función temporal externa.

```bash
python3 senal/generar_chirp.py
# Genera chirp_source.txt, que se copia a la carpeta DATA/ del caso
```

### 2. Generar la malla con Gmsh

La geometría se define en `three.geo`. El comando genera la malla en 2D y la
exporta en formato MSH versión 2, que es el que entiende el conversor de
SPECFEM2D.

```bash
gmsh three.geo -2 -format msh2 -o three.msh
```

- `-2` indica mallado bidimensional (la sección transversal del tronco).
- `-format msh2` fuerza la versión 2 del formato MSH; las versiones más recientes
  (4.x) no son compatibles con el conversor de SPECFEM2D.
- `-o three.msh` define el fichero de salida.

Antes de convertir, conviene verificar que los grupos físicos (materiales) están
bien definidos, ya que tras `BooleanFragments` la numeración de las superficies
puede cambiar:

```bash
python3 -c "
with open('three.msh') as f: content = f.read()
start = content.find('\$PhysicalNames')
end   = content.find('\$EndPhysicalNames')
print(content[start:end+17])
"
```

### 3. Convertir la malla al formato de SPECFEM2D

El script oficial de conversión transforma la malla de Gmsh en los ficheros que
SPECFEM2D necesita (Mesh_File, Nodes_File, Material_File y ficheros de borde).

```bash
python3 ~/specfem2d/utils/Gmsh/LibGmsh2Specfem_convert_Gmsh_to_Specfem2D_official.py \
        three.msh -t A -b A -l A -r A
```

- `-t A -b A -l A -r A` marcan como absorbentes (Absorbing) los cuatro bordes del
  dominio: top, bottom, left, right. Esto es necesario para que las capas PML se
  apliquen posteriormente en todo el perímetro.

Tras la conversión, se comprueban los ficheros generados:

```bash
ls -la Mesh_File Nodes_File Material_File 2>/dev/null
```

### 4. Marcar las capas absorbentes PML con convert_cpml

Las capas PML **no se definen en Gmsh**: se marcan a posteriori con la herramienta
`convert_cpml` de SPECFEM2D, que identifica los elementos del borde del dominio y
los etiqueta como PML. Por eso en la malla de Gmsh los elementos del borde son
cuadriláteros idénticos al resto del dominio.

```bash
# Enlazar los ficheros de malla con los nombres que espera convert_cpml
ln -sf Mesh_File  mesh_file
ln -sf Nodes_File nodes_coords_file

# Ejecutar convert_cpml
~/specfem2d/convert_cpml
```

Cuando el programa lo solicite, se introducen estos valores por teclado, en este
orden:

```
1        → definir el espesor de las PML manualmente
2        → aplicar PML en todos los bordes (incluido el superior)
0.051    → espesor de la PML en el borde Xmin (mínimo en X)
0.051    → espesor de la PML en el borde Xmax (máximo en X)
0.051    → espesor de la PML en el borde Zmin (mínimo en Z)
0.051    → espesor de la PML en el borde Zmax (máximo en Z)
```

**Por qué estos valores:**
- El **1** indica que se especificará el espesor de las capas de forma manual, en
  lugar de dejar que el programa lo estime automáticamente. El control manual
  garantiza que el espesor coincida con el previsto en la geometría (5 cm de PML
  sobre un dominio físico de 1,0 m, para un dominio total de 1,10 m).
- El **2** aplica las capas absorbentes en **los cuatro bordes** del dominio,
  incluido el superior. Esto es necesario porque el frente de onda esférico se
  propaga en todas las direcciones y debe absorberse en todo el perímetro para
  evitar reflexiones espurias que contaminarían el sismograma.
- Los **cuatro valores 0.051** (uno por cada borde: Xmin, Xmax, Zmin, Zmax) fijan
  el espesor de la capa absorbente en 0,051 m, ligeramente superior a los 0,05 m
  nominales para asegurar que la capa cubre al menos una fila completa de
  elementos (de tamaño 0,01 m) con un pequeño margen. Un espesor insuficiente
  dejaría reflexiones residuales; uno excesivo desperdiciaría dominio útil.

El programa genera los ficheros `absorbing_cpml_file` y `absorbing_surface_file`,
que SPECFEM2D lee durante la simulación.

### 5. Activar las PML en el Par_file

Finalmente, en el `Par_file` se activan las condiciones PML y el número de
elementos de espesor:

```bash
PML_BOUNDARY_CONDITIONS = .true.
NELEM_PML_THICKNESS     = 5
```

### 6. Ejecutar la simulación

```bash
./run_this_example.sh
```

Este script ejecuta el mallador (si procede) y el solver, y deja los resultados en
`OUTPUT_FILES/`. Cada simulación de 600 µs de propagación tarda aproximadamente
24–25 segundos.

Para los casos con múltiples receptores, se verifica que se han generado los 16
sismogramas:

```bash
ls OUTPUT_FILES/AA.S00*.BXX.semd | wc -l   # debe devolver 16
```

### 7. Postprocesar en MATLAB

Los sismogramas (`.semd`) se procesan con los scripts de MATLAB, que calculan el
tiempo de vuelo por tres métodos (correlación cruzada, umbral y envolvente de
energía), el retardo relativo (dTOF) y generan las figuras de análisis.

```matlab
run('scripts_matlab/validacion_par_file_v2.m')     % validación física base
run('scripts_matlab/validacion_malla_circular.m')  % geometría circular
run('scripts_matlab/analisis_16_sensores.m')       % localización angular
```

---

## Configuración de la fuente y los receptores

**Fuente** (fichero `SOURCE`): fuerza puntual con señal chirp externa
(`time_function_type = 8`). En la geometría circular es imprescindible orientar la
fuerza radialmente hacia el centro del tronco:

```
anglesource = 90.0    # fuerza radial (geometría circular)
```

En el punto tangente de la superficie curva, una fuerza vertical (`anglesource=0`)
excitaría el modo tangencial lento; con `anglesource=90` se excita el modo radial
rápido, que es el de interés. En la geometría rectangular, el borde plano es
tolerante a la orientación.

**Receptores** (fichero `STATIONS`): en los casos multisensor, 16 receptores
distribuidos cada 22,5° sobre una circunferencia de radio 0,145 m centrada en
(0,55, 0,55). Para que SPECFEM2D use este fichero se activa en el `Par_file`:

```
use_existing_STATIONS = .true.
```

---

## Referencia

Espinosa, L. et al. (2019). *Ultrasonic tomography for the sustainable and
rational management of standing trees*. Ultrasonics, 91, 242–251.

## Licencia

MIT License
