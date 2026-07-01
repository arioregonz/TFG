#!/bin/bash
# =========================
# Generate SPECFEM2D chirp source file
# Sensor R3α — Espinosa 2019 (Tabla 5)
# =========================

echo "Generating chirp_source.txt ..."
python3 << 'EOF'
import numpy as np

# ── Parámetros de simulación ──────────────
dt     = 2.0e-7       # debe coincidir con DT del Par_file
NSTEP  = 3000        # debe coincidir con NSTEP del Par_file

# ── Parámetros chirp R3α (Tabla 5 tesis) ─
Ts     = 45e-6        # duración del chirp: 45 µs
f0     = 22e3         # frecuencia inicial: 22 kHz
f1     = 50e3         # frecuencia final:   50 kHz
# Fc = 36 kHz, ΔF = 28 kHz

# ── Señal completa (relleno de ceros) ─────
t_full  = np.arange(0, NSTEP * dt, dt)
signal  = np.zeros(NSTEP)

# ── Chirp con ventana Gaussiana ───────────
t_chirp = np.arange(0, Ts, dt)
mu      = Ts / 2
sigma   = Ts / 6

# f(t) = (f1-f0)*t/Ts + f0  [fórmula exacta tesis]
f_t     = (f1 - f0) * t_chirp / Ts + f0
phase   = 2 * np.pi * f_t * t_chirp

# c(t) = cos(2π·f(t)·t) · Gaussiana  [tesis usa coseno]
window  = np.exp(-((t_chirp - mu)**2) / (2 * sigma**2))
chirp   = np.cos(phase) * window

# ── Insertar chirp en señal completa ──────
signal[:len(chirp)] = chirp

# ── Guardar ───────────────────────────────
np.savetxt("chirp_source.txt", np.column_stack((t_full, signal)))

print(f"chirp_source.txt creado:")
print(f"  dt      = {dt:.1e} s")
print(f"  NSTEP   = {NSTEP}")
print(f"  T_total = {NSTEP*dt*1e3:.2f} ms")
print(f"  Chirp   = {f0/1e3:.0f}–{f1/1e3:.0f} kHz, Ts={Ts*1e6:.0f} µs")
print(f"  Puntos chirp / total: {len(t_chirp)} / {NSTEP}")
EOF
echo "Done."
