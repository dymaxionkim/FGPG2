"""Tkinter GUI for FGPG2 gear profile generator — ttk flat style."""

from __future__ import annotations

import os
import sys
import tkinter as tk
from tkinter import filedialog, ttk
import tkinter.font as tkfont

import pandas as pd
from PIL import Image, ImageTk

from .gear import GearParams
from .plotter import generate

DEFAULT_PARAMS = GearParams(
    m=1.0, z=18, alpha=20, x=0.0, b=0.05, a=1.0, d=1.25, c=0.2, e=0.1,
    x_0=0.0, y_0=0.0, seg_circle=360, seg_involute=15, seg_edge_r=5,
    seg_root_r=5, seg_outer=5, seg_root=5, scale=0.7,
    gen_ratio=3.0, profile="involute",
)

_INT_FIELDS = {
    "z", "seg_circle", "seg_involute", "seg_edge_r",
    "seg_root_r", "seg_outer", "seg_root",
}

GEAR_FIELDS = [
    ("m",      "Module, m = ",           "[mm] > 0"),
    ("z",      "Teeth Number, z = ",     "[ea] : - for Internal"),
    ("alpha",  "Pressure Angle, alpha = ","[deg] : 20 for Standard"),
    ("x",      "Offset Factor, x = ",    "-1.0 ~ +1.0"),
    ("b",      "Backlash Factor, b = ",  "0.0 ~ +1.0"),
    ("a",      "Addendum Factor, a = ",  "1.0 for Standard"),
    ("d",      "Dedendum Factor, d = ",  "1.25 for Standard"),
    ("c",      "Radius of Hob end, c = ","[mm]"),
    ("e",      "Radius of Tooth end, e = ","[mm]"),
    ("gen_ratio", "Gen Ratio, = ",  "[Pd/m] : Cycloid only"),
]

GRAPHICS_FIELDS = [
    ("x_0",          "x0 = ",           "[mm]"),
    ("y_0",          "y0 = ",           "[mm]"),
    ("seg_circle",   "seg_circle = ",   "[ea]"),
    ("seg_involute", "seg_involute = ", "[ea]"),
    ("seg_edge_r",   "seg_edge_r = ",   "[ea]"),
    ("seg_root_r",   "seg_root_r = ",   "[ea]"),
    ("seg_outer",    "seg_outer = ",    "[ea]"),
    ("seg_root",     "seg_root = ",     "[ea]"),
    ("scale",        "scale = ",        "0.1 ~ 1.0"),
]

PROFILE_OPTIONS = ["Involute", "Cycloid"]
PROFILE_TO_VALUE = {"Involute": "involute", "Cycloid": "cycloid"}
VALUE_TO_PROFILE = {v: k for k, v in PROFILE_TO_VALUE.items()}


class FGPG2App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("FGPG3 — Fine Involute Gear Profile Generator")
        self.resizable(True, True)
        self.minsize(900, 600)

        if sys.platform.startswith("win"):
            self.iconbitmap("FGPG2.ico")

        self.entries: dict[str, ttk.Entry] = {}
        self.params: GearParams = DEFAULT_PARAMS
        self.work_dir: str = "./Result"

        self._configure_theme()
        self._build_ui()
        self._init_entries()

    # ------------------------------------------------------------------ theme

    def _configure_theme(self) -> None:
        self._style = ttk.Style()
        available = self._style.theme_names()
        if "vista" in available:
            self._style.theme_use("vista")
        elif "clam" in available:
            self._style.theme_use("clam")

        self._header_font = tkfont.Font(size=13, weight="bold")
        self._style.configure("Header.TLabel", font=self._header_font)

    # ------------------------------------------------------------------ UI building

    def _field_row(self, parent: ttk.Frame, key: str, label: str, hint: str, width: int = 12) -> None:
        fr = ttk.Frame(parent)
        fr.pack(fill="x", pady=1)
        ttk.Label(fr, text=label, anchor="e", width=24).pack(side="left", padx=(0, 4))
        ent = ttk.Entry(fr, width=width)
        ent.pack(side="left", padx=4)
        ttk.Label(fr, text=hint, anchor="w", width=24).pack(side="left", padx=(4, 0))
        self.entries[key] = ent

    def _build_ui(self) -> None:
        # === main container ===
        main = ttk.Frame(self, padding=(12, 10, 12, 10))
        main.pack(fill="both", expand=True)

        # === top bar ===
        top = ttk.Frame(main)
        top.pack(fill="x", pady=(0, 10))
        ttk.Label(top, text="Working Directory", anchor="e", width=20).pack(side="left", padx=(0, 6))
        self.entry_wd = ttk.Entry(top, width=48)
        self.entry_wd.pack(side="left")
        ttk.Button(top, text="Browse…", command=self._browse_wd, width=9).pack(side="left", padx=(6, 0))
        ttk.Separator(main, orient="horizontal").pack(fill="x", pady=(4, 8))

        # === content area ===
        content = ttk.Frame(main)
        content.pack(fill="both", expand=True)

        # ---- left panel ----
        left = ttk.Frame(content)
        left.pack(side="left", fill="y", padx=(0, 12))

        ttk.Label(left, text="Gear Parameters", style="Header.TLabel").pack(anchor="w")
        gear_frame = ttk.Frame(left, padding=(0, 4, 0, 0))
        gear_frame.pack(fill="x")
        for key, label, hint in GEAR_FIELDS:
            self._field_row(gear_frame, key, label, hint)

        ttk.Separator(left, orient="horizontal").pack(fill="x", pady=(12, 8))

        ttk.Label(left, text="Graphics", style="Header.TLabel").pack(anchor="w")
        graph_frame = ttk.Frame(left, padding=(0, 4, 0, 0))
        graph_frame.pack(fill="x")
        for key, label, hint in GRAPHICS_FIELDS:
            self._field_row(graph_frame, key, label, hint)

        # ---- right panel ----
        right = ttk.Frame(content)
        right.pack(side="left", fill="both", expand=True)

        self._load_default_image()
        img_container = ttk.Frame(right, padding=(0, 10, 0, 10))
        img_container.pack(fill="both", expand=True)
        self.label_image = ttk.Label(img_container, image=self._img_result)
        self.label_image.image = self._img_result
        self.label_image.pack()

        ttk.Separator(right, orient="horizontal").pack(fill="x", pady=(4, 8))

        self.status_text = ttk.Label(right, text="Ready.", anchor="w")
        self.status_text.pack(fill="x", pady=(0, 6))

        btn_frame = ttk.Frame(right)
        btn_frame.pack(fill="x")
        ttk.Button(btn_frame, text="Load", command=self._load_inputs, width=9).pack(side="left", padx=(0, 4))
        self._profile_var = tk.StringVar(value="Involute")
        ttk.Combobox(
            btn_frame, textvariable=self._profile_var,
            values=PROFILE_OPTIONS, state="readonly", width=9,
        ).pack(side="left", padx=(0, 4))
        ttk.Button(btn_frame, text="▶  Run", command=self._run, width=9).pack(side="left", padx=(0, 4))
        self._toggle_var = tk.StringVar(value="Result1")
        ttk.Checkbutton(
            btn_frame, text="Toggle image",
            command=self._toggle_image,
            variable=self._toggle_var,
            onvalue="Result1", offvalue="Result2",
        ).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="Exit", command=self.destroy, width=9).pack(side="right")

    def _load_default_image(self) -> None:
        img = Image.open("./FGPG2.png").resize((500, 500))
        self._img_result = ImageTk.PhotoImage(img)

    # ------------------------------------------------------------------ parameters

    def _init_entries(self) -> None:
        self.params = DEFAULT_PARAMS
        self.entry_wd.delete(0, "end")
        self.entry_wd.insert(0, self.work_dir)
        self._params_to_entries()

    def _params_to_entries(self) -> None:
        for key, ent in self.entries.items():
            val = getattr(self.params, key)
            ent.delete(0, "end")
            ent.insert(0, val)
        self._profile_var.set(VALUE_TO_PROFILE.get(self.params.profile, "Involute"))

    def _entries_to_params(self) -> GearParams:
        kwargs = {}
        for key, ent in self.entries.items():
            raw = ent.get()
            kwargs[key] = int(raw) if key in _INT_FIELDS else float(raw)
        kwargs["profile"] = PROFILE_TO_VALUE.get(self._profile_var.get(), "involute")
        return GearParams(**kwargs)

    def _read_wd(self) -> str:
        return self.entry_wd.get().strip()

    # ------------------------------------------------------------------ CSV I/O

    def _save_inputs_csv(self) -> None:
        p = self._entries_to_params()
        rows = [
            (key, getattr(p, key))
            for key in DEFAULT_PARAMS.__dataclass_fields__
        ]
        df = pd.DataFrame(rows, columns=["parameter", "value"])
        path = os.path.join(self.work_dir, "Inputs.csv")
        df.to_csv(path, index=False)

    def _load_inputs_csv(self) -> GearParams | None:
        path = os.path.join(self.work_dir, "Inputs.csv")
        if not os.path.exists(path):
            self.status_text.config(text=f"File not found: {path}")
            return None
        df = pd.read_csv(path, index_col="parameter")
        kwargs = {}
        for key in DEFAULT_PARAMS.__dataclass_fields__:
            if key not in df.index:
                kwargs[key] = getattr(DEFAULT_PARAMS, key)
                continue
            raw = df.loc[key, "value"]
            if key == "profile":
                kwargs[key] = str(raw).strip()
            else:
                kwargs[key] = int(raw) if key in _INT_FIELDS else float(raw)
        self.status_text.config(text=f"Loaded: {path}")
        return GearParams(**kwargs)

    # ------------------------------------------------------------------ callbacks

    def _browse_wd(self) -> None:
        wd = filedialog.askdirectory()
        if not wd:
            return
        self.entry_wd.delete(0, "end")
        self.entry_wd.insert(0, wd)
        self.work_dir = wd
        self.status_text.config(text=f"Working Directory: {wd}")

    def _load_inputs(self) -> None:
        self.work_dir = self._read_wd()
        result = self._load_inputs_csv()
        if result is not None:
            self.params = result
            self._params_to_entries()

    def _run(self) -> None:
        self.params = self._entries_to_params()
        self.work_dir = self._read_wd()
        os.makedirs(self.work_dir, exist_ok=True)
        self.status_text.config(text="Generating…")
        self.update_idletasks()
        try:
            generate(self.params, self.work_dir)
        except Exception as exc:
            self.status_text.config(text=f"Error: {exc}")
            return
        self._save_inputs_csv()
        self._display_image("Result1.png")
        self.status_text.config(text="Finished.")

    def _toggle_image(self) -> None:
        name = "Result1.png" if self._toggle_var.get() == "Result1" else "Result2.png"
        self._display_image(name)

    def _display_image(self, name: str) -> None:
        path = os.path.join(self.work_dir, name)
        if not os.path.exists(path):
            self.status_text.config(text=f"Image not found: {path}")
            return
        img = Image.open(path).resize((500, 500))
        photo = ImageTk.PhotoImage(img)
        self.label_image.configure(image=photo)
        self.label_image.image = photo
        self.status_text.config(text=f"Loaded: {name}")