# -*- coding: utf-8 -*-
"""
curso_funciones_geogebra_pro.py

Descripción:
Aplicación de escritorio para un curso interactivo de funciones.
Permite a los usuarios explorar funciones lineales y cuadráticas
con un graficador similar a GeoGebra, guardando su progreso en una
base de datos SQLite.

Mejoras incluidas:
 - Corrección del error 'AttributeError: property 'style' of 'App' object has no setter'
 - Código más legible y modularizado dentro de la clase App.
 - Manejo de excepciones más específico para errores de base de datos y archivos.
 - Mejora visual con ttkbootstrap (con fallback a ttk).

Ejecutar: python curso_funciones_geogebra_pro_mejorado.py
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import sqlite3
from datetime import datetime
import numpy as np
import os
import sys

# --- Visual (ttkbootstrap opcional) ---
USE_BOOTSTRAP = False
try:
    import ttkbootstrap as tb
    USE_BOOTSTRAP = True
except ImportError:
    # Si ttkbootstrap no está instalado, usamos ttk estándar.
    USE_BOOTSTRAP = False
except Exception as e:
    # Manejo de otros posibles errores al importar
    print(f"Advertencia: No se pudo cargar ttkbootstrap. Usando ttk estándar. Error: {e}", file=sys.stderr)
    USE_BOOTSTRAP = False

# --- Matplotlib embebido ---
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import Slider
from matplotlib.backends.backend_pdf import PdfPages

DB_FILE = "curso_funciones.db"

# -------------------------
# Base de datos (SQLite)
# -------------------------
def init_db():
    """Inicializa las tablas de la base de datos si no existen."""
    try:
        conn = sqlite3.connect(DB_FILE)
        c = conn.cursor()
        c.execute("""
        CREATE TABLE IF NOT EXISTS users(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            username TEXT UNIQUE NOT NULL,
            created_at TEXT NOT NULL
        )""")
        c.execute("""
        CREATE TABLE IF NOT EXISTS progreso(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            fecha TEXT NOT NULL,
            user_id INTEGER NOT NULL,
            unidad TEXT NOT NULL,
            completado INTEGER NOT NULL,
            puntuacion INTEGER,
            FOREIGN KEY(user_id) REFERENCES users(id)
        )""")
        conn.commit()
    except sqlite3.Error as e:
        print(f"Error de base de datos al inicializar: {e}", file=sys.stderr)
        messagebox.showerror("Error de BD", "No se pudo inicializar la base de datos. Por favor, reinicie la aplicación.")
    finally:
        conn.close()

def get_or_create_user(username: str) -> int:
    """Obtiene el ID de un usuario existente o crea uno nuevo."""
    try:
        conn = sqlite3.connect(DB_FILE)
        c = conn.cursor()
        c.execute("SELECT id FROM users WHERE username = ?", (username,))
        row = c.fetchone()
        if row:
            uid = row[0]
        else:
            c.execute("INSERT INTO users (username, created_at) VALUES (?,?)",
                      (username, datetime.now().isoformat(sep=" ", timespec="seconds")))
            conn.commit()
            uid = c.lastrowid
        return uid
    except sqlite3.Error as e:
        print(f"Error de base de datos al buscar/crear usuario: {e}", file=sys.stderr)
        messagebox.showerror("Error de BD", "No se pudo gestionar el usuario.")
        return -1 # Retorna un valor que indica un error
    finally:
        conn.close()

def save_progress(user_id: int, unidad: str, completado: bool, puntuacion=None):
    """Guarda el progreso de una unidad para un usuario."""
    if user_id == -1: # No guardar si hubo un error de login
        return
    try:
        conn = sqlite3.connect(DB_FILE)
        c = conn.cursor()
        c.execute(
            "INSERT INTO progreso (fecha, user_id, unidad, completado, puntuacion) VALUES (?,?,?,?,?)",
            (datetime.now().isoformat(sep=" ", timespec="seconds"), user_id, unidad, int(bool(completado)), puntuacion)
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"Error de base de datos al guardar progreso: {e}", file=sys.stderr)
        messagebox.showerror("Error de BD", "No se pudo guardar el progreso.")
    finally:
        conn.close()

def get_stats(user_id: int):
    """Obtiene las estadísticas de progreso para un usuario."""
    if user_id == -1:
        return []
    try:
        conn = sqlite3.connect(DB_FILE)
        c = conn.cursor()
        c.execute("SELECT unidad, COUNT(*), AVG(puntuacion) FROM progreso WHERE user_id=? GROUP BY unidad", (user_id,))
        rows = c.fetchall()
        return rows
    except sqlite3.Error as e:
        print(f"Error de base de datos al obtener estadísticas: {e}", file=sys.stderr)
        messagebox.showerror("Error de BD", "No se pudieron cargar las estadísticas.")
        return []
    finally:
        conn.close()

# -------------------------
# Contenido y preguntas
# -------------------------
LESSONS = [
    {
        "id": 1,
        "title": "Introducción a Funciones",
        "content": (
            "Una función relaciona cada elemento del dominio con exactamente un elemento del rango.\n\n"
            "Notación: $f(x) = y$\n\n"
            "Conceptos clave:\n"
            "- Dominio: valores de entrada\n"
            "- Rango: valores de salida\n"
            "- Variable independiente: $x$\n"
            "- Variable dependiente: $y$\n"
        )
    },
    {
        "id": 2,
        "title": "Funciones Lineales",
        "content": (
            "Forma general: $f(x) = mx + b$\n\n"
            "- $m$: pendiente\n"
            "- $b$: ordenada al origen ($f(0)$)\n\n"
            "Características:\n"
            "- Gráfica: recta\n"
            "- Tasa de cambio constante\n\n"
            "Ejemplo: $f(x) = 3x + 2$ -> pendiente 3, ordenada 2\n"
        )
    },
    {
        "id": 3,
        "title": "Funciones Cuadráticas",
        "content": (
            "Forma general: $f(x) = ax^2 + bx + c$\n\n"
            "- $a$ determina apertura ($a>0$: arriba, $a<0$: abajo)\n"
            "- Vértice: $x_v = -b/(2a)$, $y_v = f(x_v)$\n"
            "- Raíces: soluciones de $f(x)=0$ (fórmula cuadrática o factorización)\n"
            "\nEjemplo: $f(x) = x^2 - 4x + 3$"
        )
    },
    {
        "id": 4,
        "title": "Evaluación Final",
        "content": "Cuestionario final para evaluar lo aprendido (se guarda la puntuación)."
    }
]

QUIZ = [
    {"q": "¿Cuál es la forma general de una función lineal?", "opts": ["$y = ax^2 + bx + c$", "$y = mx + b$", "$y = x$"], "a": "$y = mx + b$"},
    {"q": "En $f(x)=3x+5$, ¿cuál es la pendiente?", "opts": ["5", "3", "8"], "a": "3"},
    {"q": "La gráfica de una función cuadrática es:", "opts": ["Recta", "Círculo", "Parábola"], "a": "Parábola"},
    {"q": "Si $f(x)=2x^2-8x+6$, la abscisa del vértice es:", "opts": ["2", "-2", "4"], "a": "2"}
]

# -------------------------
# Graficadora interactiva
# -------------------------
class InteractiveGrapher:
    """
    Clase que gestiona la gráfica interactiva (Matplotlib embebido).
    Permite manipular funciones lineales y cuadráticas con sliders y puntos
    arrastrables.
    """
    def __init__(self, tk_parent):
        self.tk_parent = tk_parent
        self.mode = "linear"
        self.control_mode = "sliders"
        self.params = {"m": 1.0, "b": 0.0, "a": 1.0, "bb": 0.0, "c": 0.0}

        self.fig = plt.Figure(figsize=(8, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self._style_axes()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.tk_parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.tk_parent, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X, pady=2)

        self.slider_axes = []
        self.sliders = {}
        self._create_slider_axes()
        self._build_sliders_for_mode("linear")

        self.line_plot = None
        self.scatter_controls = []
        self.dragging_idx = None

        self._connect_events()
        self._update_plot()

    def _style_axes(self):
        self.ax.clear()
        self.ax.axhline(0, color="#222", linewidth=0.8)
        self.ax.axvline(0, color="#222", linewidth=0.8)
        self.ax.set_xlim(-10, 10)
        self.ax.set_ylim(-10, 10)
        self.ax.grid(True, linestyle="--", alpha=0.5)
        self.ax.set_title("Graficadora Interactiva (tipo GeoGebra)")

    def set_mode(self, mode: str):
        """Cambia el modo de la función a graficar (lineal/cuadrática)."""
        assert mode in ("linear", "quadratic")
        self.mode = mode
        self._style_axes()
        self._clear_control_points()
        self._rebuild_sliders_for_mode()
        self._update_plot()

    def set_control_mode(self, how: str):
        """Cambia el método de control (sliders/puntos arrastrables)."""
        assert how in ("sliders", "puntos")
        self.control_mode = how
        self._style_axes()
        self._clear_control_points()
        if how == "puntos":
            if self.mode == "linear":
                pts = [(-3, self.params["m"] * (-3) + self.params["b"]),
                       (3, self.params["m"] * (3) + self.params["b"])]
            else:
                def f(x):
                    return self.params["a"] * x**2 + self.params["bb"] * x + self.params["c"]
                pts = [(-3, f(-3)), (0, f(0)), (3, f(3))]
            for (x, y) in pts:
                sc = self.ax.scatter([x], [y], s=80, color="#ff7f0e", zorder=3, picker=5)
                self.scatter_controls.append(sc)
        self._update_plot()

    # ---------- Sliders ----------
    def _create_slider_axes(self):
        """Crea los ejes de los sliders en la parte inferior de la figura."""
        for sa in self.slider_axes:
            sa.remove()
        self.slider_axes = []
        bottom = 0.01
        height = 0.03
        self.slider_axes.append(self.fig.add_axes([0.12, bottom, 0.78, height]))
        self.slider_axes.append(self.fig.add_axes([0.12, bottom + 0.04, 0.78, height]))
        self.slider_axes.append(self.fig.add_axes([0.12, bottom + 0.08, 0.78, height]))

    def _rebuild_sliders_for_mode(self):
        """Reconstruye los sliders según el modo de función actual."""
        for name, sl in list(self.sliders.items()):
            try:
                sl.ax.clear()
            except Exception:
                pass
        self.sliders = {}
        self._create_slider_axes()
        self._build_sliders_for_mode(self.mode)

    def _build_sliders_for_mode(self, mode):
        """Crea los widgets de slider para el modo especificado."""
        for ax in self.slider_axes:
            ax.clear()
        self.sliders = {}

        if mode == "linear":
            ax_m, ax_b = self.slider_axes[2], self.slider_axes[1]
            self.sliders["m"] = Slider(ax=ax_m, label="m (pendiente)", valmin=-5, valmax=5, valinit=self.params["m"], valstep=0.05)
            self.sliders["b"] = Slider(ax=ax_b, label="b (ordenada)", valmin=-10, valmax=10, valinit=self.params["b"], valstep=0.1)

            self.sliders["m"].on_changed(lambda v: self._on_slider_change("m", v))
            self.sliders["b"].on_changed(lambda v: self._on_slider_change("b", v))
        else:
            ax_a, ax_b, ax_c = self.slider_axes[2], self.slider_axes[1], self.slider_axes[0]
            self.sliders["a"] = Slider(ax=ax_a, label="a", valmin=-3, valmax=3, valinit=self.params["a"], valstep=0.01)
            self.sliders["bb"] = Slider(ax=ax_b, label="b", valmin=-10, valmax=10, valinit=self.params["bb"], valstep=0.05)
            self.sliders["c"] = Slider(ax=ax_c, label="c", valmin=-10, valmax=10, valinit=self.params["c"], valstep=0.1)

            self.sliders["a"].on_changed(lambda v: self._on_slider_change("a", v))
            self.sliders["bb"].on_changed(lambda v: self._on_slider_change("bb", v))
            self.sliders["c"].on_changed(lambda v: self._on_slider_change("c", v))

        self.fig.canvas.draw_idle()

    def _on_slider_change(self, name, value):
        """Actualiza el parámetro al mover un slider."""
        self.params[name] = float(value)
        if self.control_mode == "puntos":
            self._reposition_points_from_params()
        self._update_plot()

    # ---------- Control por puntos ----------
    def _clear_control_points(self):
        """Elimina los puntos de control arrastrables de la gráfica."""
        for sc in self.scatter_controls:
            try:
                sc.remove()
            except Exception:
                pass
        self.scatter_controls = []
        self.dragging_idx = None

    def _reposition_points_from_params(self):
        """
        Calcula las nuevas posiciones de los puntos de control en base a los parámetros
        actuales de la función.
        """
        if not self.scatter_controls:
            return
        if self.mode == "linear":
            xs = [-3, 3]
            ys = [self.params["m"] * xs[0] + self.params["b"],
                  self.params["m"] * xs[1] + self.params["b"]]
        else:
            def f(x): return self.params["a"] * x**2 + self.params["bb"] * x + self.params["c"]
            xs = [-3, 0, 3]
            ys = [f(xs[0]), f(xs[1]), f(xs[2])]
        for sc, x, y in zip(self.scatter_controls, xs, ys):
            sc.set_offsets([[x, y]])

    def _connect_events(self):
        """Conecta los eventos del mouse a las funciones de interacción."""
        self.cid_press = self.fig.canvas.mpl_connect("button_press_event", self._on_press)
        self.cid_release = self.fig.canvas.mpl_connect("button_release_event", self._on_release)
        self.cid_motion = self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion)

    def _on_press(self, event):
        """Detecta si se hizo clic en un punto de control."""
        if event.inaxes != self.ax or self.control_mode != "puntos":
            return
        for i, sc in enumerate(self.scatter_controls):
            if sc.get_paths() and sc.contains(event)[0]:
                self.dragging_idx = i
                break

    def _on_motion(self, event):
        """Mueve el punto de control y actualiza la función en tiempo real."""
        if self.dragging_idx is None or event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        sc = self.scatter_controls[self.dragging_idx]
        sc.set_offsets([[event.xdata, event.ydata]])
        self._recompute_params_from_points()
        self._sync_sliders_from_params()
        self._update_plot(light=True)

    def _on_release(self, event):
        """Finaliza el arrastre del punto de control."""
        self.dragging_idx = None
        self._update_plot()

    def _recompute_params_from_points(self):
        """
        Recalcula los parámetros de la función (m, b, a, c) en base a
        la posición de los puntos de control arrastrables.
        """
        if not self.scatter_controls:
            return
        pts = [tuple(sc.get_offsets()[0]) for sc in self.scatter_controls]
        if self.mode == "linear" and len(pts) >= 2:
            (x1, y1), (x2, y2) = pts[0], pts[1]
            if abs(x2 - x1) < 1e-8:
                return
            m = (y2 - y1) / (x2 - x1)
            b = y1 - m * x1
            self.params["m"] = float(m)
            self.params["b"] = float(b)
        elif self.mode == "quadratic" and len(pts) >= 3:
            xs = np.array([p[0] for p in pts[:3]])
            ys = np.array([p[1] for p in pts[:3]])
            A = np.vstack([xs**2, xs, np.ones_like(xs)]).T
            try:
                sol = np.linalg.solve(A, ys)
                self.params["a"], self.params["bb"], self.params["c"] = [float(s) for s in sol]
            except np.linalg.LinAlgError:
                pass

    def _sync_sliders_from_params(self):
        """Sincroniza los sliders para que coincidan con los parámetros recalculados."""
        for k, v in self.params.items():
            if k in self.sliders:
                try:
                    self.sliders[k].set_val(v)
                except Exception:
                    pass

    # ---------- Plot principal ----------
    def _update_plot(self, light=False):
        """
        Dibuja la función y los puntos de control.
        `light=True` se usa para actualizaciones de arrastre, evitando
        redibujar componentes estáticos para una mejor performance.
        """
        self.ax.collections = [c for c in self.ax.collections if c in self.scatter_controls]
        self.ax.lines = []
        self.ax.texts = []
        self._style_axes()

        x = np.linspace(-10, 10, 600)
        if self.mode == "linear":
            m, b = self.params["m"], self.params["b"]
            y = m * x + b
            label = f"$y = {m:.2f}x + {b:.2f}$"
        else:
            a, bb, c = self.params["a"], self.params["bb"], self.params["c"]
            y = a * x**2 + bb * x + c
            label = f"$y = {a:.2f}x² + {bb:.2f}x + {c:.2f}$"

        self.ax.plot(x, y, lw=2, color="#2563eb", label=label)
        self.ax.legend(loc="upper left")

        for sc in self.scatter_controls:
            self.ax.add_collection(sc)

        self.fig.canvas.draw_idle()

    # ---------- Exportaciones ----------
    def save_png(self, path):
        """Guarda la gráfica actual como una imagen PNG."""
        self.fig.savefig(path, dpi=300, bbox_inches="tight")

    def save_pdf(self, path):
        """Guarda la gráfica actual como un archivo PDF."""
        self.fig.savefig(path, dpi=300, bbox_inches="tight")

    def export_report_pdf(self, path, username: str, lesson_title: str):
        """Genera un reporte PDF con resumen y gráfica."""
        with PdfPages(path) as pdf:
            fig1 = plt.figure(figsize=(8.27, 11.69))
            t = fig1.text
            t(0.1, 0.9, "Reporte de Progreso", fontsize=18, weight="bold")
            t(0.1, 0.85, f"Usuario: {username}", fontsize=12)
            t(0.1, 0.82, f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", fontsize=12)
            t(0.1, 0.79, f"Lección: {lesson_title}", fontsize=12)
            t(0.1, 0.74, "Parámetros actuales:", fontsize=12, weight="bold")

            if self.mode == "linear":
                m, b = self.params["m"], self.params["b"]
                t(0.1, 0.70, f"Modo: Lineal — $y = {m:.4f}x + {b:.4f}$", fontsize=11)
            else:
                a, bb, c = self.params["a"], self.params["bb"], self.params["c"]
                t(0.1, 0.70, f"Modo: Cuadrática — $y = {a:.4f}x² + {bb:.4f}x + {c:.4f}$", fontsize=11)

            t(0.1, 0.64, "Notas:", fontsize=12, weight="bold")
            t(0.1, 0.61, "• Interacción tipo GeoGebra con sliders y puntos de control.", fontsize=10)
            t(0.1, 0.58, "• Exportado desde la aplicación Python (Tkinter + Matplotlib).", fontsize=10)
            pdf.savefig(fig1)
            plt.close(fig1)

            fig2 = plt.figure(figsize=(8.27, 11.69))
            ax2 = fig2.add_subplot(111)
            ax2.axhline(0, color="#222", linewidth=0.8)
            ax2.axvline(0, color="#222", linewidth=0.8)
            ax2.set_xlim(-10, 10)
            ax2.set_ylim(-10, 10)
            ax2.grid(True, linestyle="--", alpha=0.5)
            x = np.linspace(-10, 10, 600)
            if self.mode == "linear":
                m, b = self.params["m"], self.params["b"]
                y = m * x + b
                label = f"$y = {m:.2f}x + {b:.2f}$"
            else:
                a, bb, c = self.params["a"], self.params["bb"], self.params["c"]
                y = a * x**2 + bb * x + c
                label = f"$y = {a:.2f}x² + {bb:.2f}x + {c:.2f}$"
            ax2.plot(x, y, lw=2, color="#2563eb", label=label)
            for sc in self.scatter_controls:
                px, py = sc.get_offsets()[0]
                ax2.scatter([px], [py], s=80, color="#ff7f0e", zorder=3)
            ax2.legend(loc="upper left")
            pdf.savefig(fig2)
            plt.close(fig2)

# -------------------------
# App Principal
# -------------------------
class App(tk.Tk if not USE_BOOTSTRAP else tb.Window):
    def __init__(self):
        # Corregido: Se pasa el tema al constructor en lugar de asignar a `self.style`
        if USE_BOOTSTRAP:
            super().__init__(themename="litera", title="Curso: Funciones (GeoGebra-like)")
        else:
            super().__init__(title="Curso: Funciones (GeoGebra-like)")
        
        self.geometry("1280x800")
        init_db()

        self.user_id = None
        self.username = None

        self.current_index = 0
        self.completed = set()
        self.user_answers = {}

        self._login_dialog()
        # Si el login falló, salimos
        if self.user_id == -1:
            self.quit()
            return
        
        self._build_ui()
        self.load_lesson(0)

    # ---------- Login ----------
    def _login_dialog(self):
        top = tk.Toplevel(self)
        top.title("Login de Usuario")
        top.grab_set()
        
        # Simplificación: define las clases una vez
        WidgetFrame = tb.Frame if USE_BOOTSTRAP else ttk.Frame
        WidgetLabel = tb.Label if USE_BOOTSTRAP else ttk.Label
        WidgetEntry = tb.Entry if USE_BOOTSTRAP else ttk.Entry
        WidgetButton = tb.Button if USE_BOOTSTRAP else ttk.Button

        frm = WidgetFrame(top, padding=15)
        frm.pack(fill=tk.BOTH, expand=True)
        
        WidgetLabel(frm, text="Ingrese su nombre de usuario:").pack(anchor="w", pady=(0, 6))
        username_var = tk.StringVar()
        ent = WidgetEntry(frm, textvariable=username_var, width=30)
        ent.pack(anchor="w")
        ent.focus_set()

        def do_login():
            name = username_var.get().strip()
            if not name:
                messagebox.showwarning("Atención", "Ingrese un nombre de usuario.")
                return
            uid = get_or_create_user(name)
            self.user_id = uid
            self.username = name
            top.destroy()

        WidgetButton(frm, text="Ingresar", command=do_login, bootstyle="primary" if USE_BOOTSTRAP else None).pack(pady=10)
        self.wait_window(top)

    # ---------- UI ----------
    def _build_ui(self):
        # Simplificación: define las clases una vez
        WidgetFrame = tb.Frame if USE_BOOTSTRAP else ttk.Frame
        WidgetLabel = tb.Label if USE_BOOTSTRAP else ttk.Label
        WidgetButton = tb.Button if USE_BOOTSTRAP else ttk.Button
        WidgetRadiobutton = tb.Radiobutton if USE_BOOTSTRAP else ttk.Radiobutton

        main = WidgetFrame(self)
        main.pack(fill=tk.BOTH, expand=True)

        side = WidgetFrame(main, width=280)
        side.pack(side=tk.LEFT, fill=tk.Y)
        side.pack_propagate(False)

        title = WidgetLabel(side, text="📚 Contenido del Curso", font=("Inter", 14, "bold"))
        title.pack(padx=12, pady=(12, 6))

        self.lesson_buttons = []
        for i, l in enumerate(LESSONS):
            b = WidgetButton(side, text=l["title"], command=lambda i=i: self.load_lesson(i),
                             bootstyle="info-outline" if USE_BOOTSTRAP else None)
            b.pack(fill=tk.X, padx=12, pady=6)
            self.lesson_buttons.append(b)

        center = WidgetFrame(main)
        center.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.content_frame = WidgetFrame(center, padding=10)
        self.content_frame.pack(fill=tk.X, padx=10, pady=10)

        graph_container = WidgetFrame(center, padding=8)
        graph_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))

        top_controls = WidgetFrame(graph_container)
        top_controls.pack(fill=tk.X)

        self.mode_var = tk.StringVar(value="linear")
        rb1 = WidgetRadiobutton(top_controls, text="Lineal", variable=self.mode_var, value="linear", command=self._on_mode_change)
        rb2 = WidgetRadiobutton(top_controls, text="Cuadrática", variable=self.mode_var, value="quadratic", command=self._on_mode_change)
        rb1.pack(side=tk.LEFT, padx=6, pady=4)
        rb2.pack(side=tk.LEFT, padx=6, pady=4)

        self.ctrl_var = tk.StringVar(value="sliders")
        rbS = WidgetRadiobutton(top_controls, text="Sliders", variable=self.ctrl_var, value="sliders", command=self._on_ctrl_change)
        rbP = WidgetRadiobutton(top_controls, text="Puntos", variable=self.ctrl_var, value="puntos", command=self._on_ctrl_change)
        rbS.pack(side=tk.LEFT, padx=6)
        rbP.pack(side=tk.LEFT, padx=6)

        btn_png = WidgetButton(top_controls, text="Exportar PNG", bootstyle="success" if USE_BOOTSTRAP else None, command=self._export_png)
        btn_pdf = WidgetButton(top_controls, text="Exportar PDF", bootstyle="primary" if USE_BOOTSTRAP else None, command=self._export_pdf)
        btn_rep = WidgetButton(top_controls, text="Reporte PDF", bootstyle="warning" if USE_BOOTSTRAP else None, command=self._export_report)
        btn_png.pack(side=tk.RIGHT, padx=4)
        btn_pdf.pack(side=tk.RIGHT, padx=4)
        btn_rep.pack(side=tk.RIGHT, padx=4)

        self.grapher = InteractiveGrapher(graph_container)

        right = WidgetFrame(main, width=280)
        right.pack(side=tk.RIGHT, fill=tk.Y)
        right.pack_propagate(False)

        lbl = WidgetLabel(right, text="Estadísticas", font=("Inter", 12, "bold"))
        lbl.pack(pady=10)

        self.stats_text = tk.Text(right, height=25, width=30)
        self.stats_text.pack(padx=8)
        self._refresh_stats()

    # ---------- Carga de lecciones ----------
    def load_lesson(self, idx):
        self.current_index = idx
        lesson = LESSONS[idx]
        for w in self.content_frame.winfo_children():
            w.destroy()

        WidgetLabel = tb.Label if USE_BOOTSTRAP else ttk.Label
        WidgetFrame = tb.Frame if USE_BOOTSTRAP else ttk.Frame
        WidgetButton = tb.Button if USE_BOOTSTRAP else ttk.Button
        WidgetRadiobutton = tb.Radiobutton if USE_BOOTSTRAP else ttk.Radiobutton

        title = WidgetLabel(self.content_frame, text=lesson["title"], font=("Inter", 16, "bold"))
        title.pack(anchor="w")

        txt = tk.Text(self.content_frame, wrap="word", height=10, bd=0)
        txt.insert("1.0", lesson["content"])
        txt.config(state="disabled")
        txt.pack(fill=tk.X, pady=(6, 10))

        act = WidgetFrame(self.content_frame)
        act.pack(fill=tk.X, pady=4)

        if lesson["id"] == 1:
            self._build_intro_activity(act)
        elif lesson["id"] == 2:
            self._build_linear_activity(act)
            self.mode_var.set("linear")
            self._on_mode_change()
            self.ctrl_var.set("sliders")
            self._on_ctrl_change()
        elif lesson["id"] == 3:
            self._build_quadratic_activity(act)
            self.mode_var.set("quadratic")
            self._on_mode_change()
            self.ctrl_var.set("sliders")
            self._on_ctrl_change()
        else:
            self._build_quiz_activity(act)

        nav = WidgetFrame(self.content_frame)
        nav.pack(fill=tk.X, pady=6)
        if idx > 0:
            WidgetButton(nav, text="← Anterior", command=lambda: self.load_lesson(idx - 1)).pack(side=tk.LEFT, padx=4)
        WidgetButton(nav, text="Completar ✓", command=self._mark_completed).pack(side=tk.LEFT, padx=4)
        if idx < len(LESSONS) - 1:
            WidgetButton(nav, text="Siguiente →", command=lambda: self.load_lesson(idx + 1)).pack(side=tk.LEFT, padx=4)

        self._refresh_stats()

    def _build_intro_activity(self, parent):
        WidgetLabel = tb.Label if USE_BOOTSTRAP else ttk.Label
        lbl = WidgetLabel(parent, text="Actividad: Escribe con tus palabras qué es una función.")
        lbl.pack(anchor="w")
        self.reflection = tk.Text(parent, height=4)
        self.reflection.pack(fill=tk.X, pady=4)

    def _build_linear_activity(self, parent):
        WidgetLabel = tb.Label if USE_BOOTSTRAP else ttk.Label
        WidgetLabel(parent, text="Explora la recta con sliders o arrastrando 2 puntos en la gráfica.").pack(anchor="w")

    def _build_quadratic_activity(self, parent):
        WidgetLabel = tb.Label if USE_BOOTSTRAP else ttk.Label
        WidgetLabel(parent, text="Explora la parábola con sliders o arrastrando 3 puntos en la gráfica.").pack(anchor="w")

    def _build_quiz_activity(self, parent):
        WidgetLabel = tb.Label if USE_BOOTSTRAP else ttk.Label
        WidgetRadiobutton = tb.Radiobutton if USE_BOOTSTRAP else ttk.Radiobutton
        WidgetButton = tb.Button if USE_BOOTSTRAP else ttk.Button

        WidgetLabel(parent, text="Responde y presiona Enviar. Se guarda la puntuación.").pack(anchor="w")
        self.q_vars = []
        for i, q in enumerate(QUIZ):
            WidgetLabel(parent, text=f"{i + 1}. {q['q']}").pack(anchor="w", pady=(8, 2))
            var = tk.StringVar(value="")
            self.q_vars.append(var)
            for opt in q['opts']:
                WidgetRadiobutton(parent, text=opt, variable=var, value=opt).pack(anchor="w")
        WidgetButton(parent, text="Enviar", bootstyle="success" if USE_BOOTSTRAP else None, command=self._submit_quiz).pack(pady=8)

    # ---------- Nav / progreso / stats ----------
    def _mark_completed(self):
        lesson = LESSONS[self.current_index]
        self.completed.add(lesson["id"])
        save_progress(self.user_id, lesson["title"], True, None)
        messagebox.showinfo("Completado", f"Lección completada: {lesson['title']}")
        self._refresh_stats()

    def _submit_quiz(self):
        score = sum(1 for i, q in enumerate(QUIZ) if self.q_vars[i].get() == q['a'])
        save_progress(self.user_id, "Cuestionario Final", True, score)
        messagebox.showinfo("Resultado", f"Puntuación: {score}/{len(QUIZ)}")
        self._refresh_stats()

    def _refresh_stats(self):
        rows = get_stats(self.user_id)
        self.stats_text.config(state="normal")
        self.stats_text.delete("1.0", tk.END)
        self.stats_text.insert(tk.END, f"Usuario: {self.username}\n\n")
        if not rows:
            self.stats_text.insert(tk.END, "Sin datos aún.\n")
        else:
            for unidad, count, avg in rows:
                puntuacion_str = f"{avg:.2f}" if avg is not None else '-'
                self.stats_text.insert(tk.END, f"{unidad}\n  Registros: {count}\n  Puntuación media: {puntuacion_str}\n\n")
        self.stats_text.config(state="disabled")

    # ---------- Cambios de modo/control ----------
    def _on_mode_change(self):
        self.grapher.set_mode(self.mode_var.get())

    def _on_ctrl_change(self):
        self.grapher.set_control_mode(self.ctrl_var.get())

    # ---------- Exportaciones ----------
    def _export_png(self):
        f = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("Imagen PNG", "*.png")],
            initialfile="grafica.png"
        )
        if not f: return
        try:
            self.grapher.save_png(f)
            messagebox.showinfo("Exportar PNG", f"Imagen guardada en:\n{f}")
        except Exception as e:
            messagebox.showerror("Error", f"Error al guardar el archivo: {e}")

    def _export_pdf(self):
        f = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF", "*.pdf")],
            initialfile="grafica.pdf"
        )
        if not f: return
        try:
            self.grapher.save_pdf(f)
            messagebox.showinfo("Exportar PDF", f"PDF guardado en:\n{f}")
        except Exception as e:
            messagebox.showerror("Error", f"Error al guardar el archivo: {e}")

    def _export_report(self):
        f = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF", "*.pdf")],
            initialfile="reporte.pdf"
        )
        if not f: return
        try:
            lesson_title = LESSONS[self.current_index]["title"]
            self.grapher.export_report_pdf(f, self.username, lesson_title)
            messagebox.showinfo("Reporte", f"Reporte PDF guardado en:\n{f}")
        except Exception as e:
            messagebox.showerror("Error", f"Error al guardar el reporte: {e}")

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    app = App()
    app.mainloop()
