# Caracterización de ondas químicas en flujo laminar (BZ - Poiseuille)

---
![Flujo laminar de Poiseuille](figs/poiseulle.jpg)
## 🧪 Descripción general

Este proyecto tiene como objetivo estudiar el comportamiento de las ondas químicas generadas por la reacción de Belousov-Zhabotinsky (BZ) en presencia de un **flujo laminar tipo Poiseuille** dentro de una **caja rectangular** con:

- Condiciones **periódicas** en los bordes izquierdo y derecho.
- Condiciones de **no flujo** en la tapa y la base.

![Esquema de condiciones de frontera](figs/c_frontera.jpg)


El modelo se basa en el **Oregonator de dos variables**, incluyendo términos de advección y una implementación numérica en **Fortran**, con análisis posterior en **Python**.

---

## 🎯 Objetivos específicos

- Analizar la **velocidad de propagación**, deformación y estabilidad de los frentes químicos bajo flujo laminar.
- Comparar diferentes tamaños de caja y condiciones iniciales para identificar simetrías, inestabilidades y patrones emergentes.
- Implementar condiciones iniciales coherentes con pulsos químicos autocatalíticos en un entorno numérico estable.

---

## 🔧 Herramientas utilizadas

- Lenguaje de simulación: **Fortran 90** (modelo Oregonator + advección)
- Análisis y visualización: **Python 3**, con bibliotecas como `numpy`, `matplotlib`, `h5py`
- Formato de salida de datos: **HDF5**
- Documentación: **Markdown, LaTeX**
- Control de versiones: **Git + GitHub**

---

## ✅ Checklist de ejecución (Fortran + HDF5 + OpenMP + Python)

### 0) Preparación (una vez)
- [ ] Instalar dependencias del sistema: `gfortran`, `libhdf5-dev`, `hdf5-tools`, `ffmpeg`
- [ ] Crear/activar entorno Python y librerías: `numpy`, `matplotlib`, `h5py`

---

### 1) Preparar caso por malla (N)
- [ ] Copiar `Poiseuille.f90` al folder del caso (ej. `src/N=400/`)
- [ ] Editar `module parameters` y fijar:
  - [ ] `N = ...`
  - [ ] `M = ...`
  - [ ] `Nsnap` y `Nsnap_total` según el experimento
  - [ ] rango de `vl` (`vl_initial`, `vl_max`, `vl_step`)
- [ ] Verificar que `io_hdf5.f90` esté en la misma carpeta del caso

---

### 2) Compilar (HDF5 + OpenMP)
Dentro de `src/N=XXX/`:
- [ ] Compilar:
  - [ ] `h5fc -O3 -fopenmp Poiseuille.f90 io_hdf5.f90 -o bz`
- [ ] Confirmar ejecutable:
  - [ ] `ls -lh bz`

---

### 3) Ejecutar simulación (OpenMP)
Dentro de `src/N=XXX/`:
- [ ] Correr con hilos:
  - [ ] `OMP_NUM_THREADS=8 ./bz`
- [ ] (Opcional) Afinidad recomendada:
  - [ ] `OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores ./bz`

---

### 4) Verificar que se generaron salidas por cada `vl`
Dentro de `src/N=XXX/`:
- [ ] Se crearon carpetas `vX_YY/`
- [ ] En un caso de prueba (ej. `v0_00/`) existen:
  - [ ] `run.h5`
  - [ ] `xcm_x.dat`

---

### 5) Copiar scripts de análisis (por cada caso N)
En `src/N=XXX/`:
- [ ] Copiar:
  - [ ] `verificador.py`
  - [ ] `velocidad.py`
  - [ ] `make_video_h5.py`

---

### 6) Verificación del HDF5 (rápido)
Dentro de `src/N=XXX/`:
- [ ] Revisar estructura:
  - [ ] `h5ls -r v0_00/run.h5`
- [ ] Verificar con Python:
  - [ ] `python3 verificador.py v0_00/run.h5`

---

### 7) Velocidad del pulso (desde `xcm_x.dat`)
Dentro de `src/N=XXX/vX_YY/`:
- [ ] Ejecutar:
  - [ ] `python3 ../velocidad.py`
- [ ] Confirmar:
  - [ ] gráfica `x_cm(t)` (unwrapped)
  - [ ] velocidad por bloques (ventanas)
  - [ ] valor final de velocidad (comparar con caso libre ~17.4)

---

### 8) Video desde HDF5 (heatmap + perfil 1D)
Dentro de `src/N=XXX/vX_YY/`:
- [ ] Video de `u`:
  - [ ] `python3 ../make_video_h5.py run.h5 --field u --j 11 --fps 20 --figsize 12 6`
- [ ] Video de `v`:
  - [ ] `python3 ../make_video_h5.py run.h5 --field v --j 11 --fps 20 --figsize 12 6`
- [ ] (Opcional) reducir peso:
  - [ ] `python3 ../make_video_h5.py run.h5 --field u --j 11 --stride 3 --fps 20 --figsize 12 6`
- [ ] (Opcional) GIF si no hay ffmpeg:
  - [ ] `python3 ../make_video_h5.py run.h5 --field u --gif --stride 3 --fps 10`

---

### 9) Checklist final por cada `vl`
- [ ] `run.h5` tiene datasets `/t`, `/u`, `/v`
- [ ] `xcm_x.dat` se llena durante la corrida
- [ ] La velocidad estimada no es cero (caso libre ≈ 17.4)
- [ ] El video muestra pulso estable (heatmap + perfil 1D coherente)


## 📁 Estructura del repositorio

```bash
caracterizacion-bz-poiseuille/
├── README.md           # Esta descripción general del proyecto
├── TODO.md             # Tareas pendientes, mejoras y bugs conocidos
├── diario.md           # Bitácora personal de desarrollo y hallazgos
├── CHANGELOG.md        # Historial de cambios estructurados
├── LICENSE             # Licencia del proyecto (MIT)
├── src/                # Códigos fuente en Fortran
├── scripts/            # Scripts en Python para análisis y visualización
├── output/             # Archivos de salida (datos simulados) - NO SE SUBEN a GitHub
├── figs/               # Figuras y gráficas generadas
├── doc/                # Documentos, borradores, artículos científicos
└── notebooks/          # Jupyter Notebooks explicativos y análisis interactivo


