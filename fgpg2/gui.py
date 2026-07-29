"""Tkinter GUI for FGPG2 gear profile generator."""

from __future__ import annotations

import os
import sys
import tkinter as tk
from tkinter import filedialog
import tkinter.font

import pandas as pd
from PIL import Image, ImageTk

from .gear import GearParams
from .plotter import generate

PADX = 1
PADY = 1

DEFAULT_PARAMS = GearParams(
    m=1.0, z=18, alpha=20, x=0.0, b=0.05, a=1.0, d=1.25, c=0.2, e=0.1,
    x_0=0.0, y_0=0.0, seg_circle=360, seg_involute=15, seg_edge_r=5,
    seg_root_r=5, seg_outer=5, seg_root=5, scale=0.7,
)

_INT_FIELDS = {"z", "seg_circle", "seg_involute", "seg_edge_r", "seg_root_r", "seg_outer", "seg_root"}

LHS_FIELDS = [
    (1,  "m",      "Module, m = ",          "[mm] > 0"),
    (2,  "z",      "Teeth Number, z = ",    "[ea] : - for Internal"),
    (3,  "alpha",  "Pressure Angle, alpha = ", "[deg] : 20 for Standard"),
    (4,  "x",      "Offset Factor, x = ",   "-1.0 ~ +1.0"),
    (5,  "b",      "Backlash Factor, b = ", "0.0 ~ +1.0"),
    (6,  "a",      "Addendum Factor, a = ", "1.0 for Standard"),
    (7,  "d",      "Dedendum Factor, d = ", "1.25 for Standard"),
    (8,  "c",      "Radius of Hob end, c = ", "[mm]"),
    (9,  "e",      "Radius of Tooth end, e = ", "[mm]"),
]

RHS_FIELDS = [
    (11, "x_0",         "x0 = ",          "[mm]"),
    (12, "y_0",         "y0 = ",          "[mm]"),
    (13, "seg_circle",  "seg_circle = ",  "[ea]"),
    (14, "seg_involute", "seg_involute = ", "[ea]"),
    (15, "seg_edge_r",  "seg_edge_r = ",  "[ea]"),
    (16, "seg_root_r",  "seg_root_r = ",  "[ea]"),
    (17, "seg_outer",   "seg_outer = ",   "[ea]"),
    (18, "seg_root",    "seg_root = ",    "[ea]"),
    (19, "scale",       "scale = ",       "0.1 ~ 1.0"),
]


class FGPG2App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("FGPG3 with tkinter")
        self.resizable(True, True)
        self.font16 = tkinter.font.Font(size=16)
        if sys.platform.startswith("win"):
            self.iconbitmap("FGPG2.ico")

        self.entries: dict[str, tk.Entry] = {}
        self.params: GearParams = DEFAULT_PARAMS
        self.work_dir: str = "./Result"

        self._build_ui()
        self._init_entries()

    # ------------------------------------------------------------------ UI building

    def _add_field(self, row: int, key: str, label_text: str, hint_text: str) -> None:
        lbl = tk.Label(self, text=label_text, anchor="e")
        lbl.grid(row=row, column=0, padx=PADX, pady=PADY, sticky="e")
        ent = tk.Entry(self)
        ent.grid(row=row, column=1, padx=PADX, pady=PADY)
        hint = tk.Label(self, text=hint_text, anchor="w")
        hint.grid(row=row, column=2, padx=PADX, pady=PADY, sticky="w")
        self.entries[key] = ent

    def _build_ui(self) -> None:
        # ---- left column (gear spec) ----
        gs = tk.Label(self, text="# Gear Spec", font=self.font16)
        gs.grid(row=0, column=0, padx=PADX, pady=PADY, sticky="w")
        for row, key, label, hint in LHS_FIELDS:
            self._add_field(row, key, label, hint)

        # ---- right column (graphics) ----
        gr = tk.Label(self, text="# Graphics", font=self.font16)
        gr.grid(row=10, column=0, padx=PADX, pady=PADY, sticky="w")
        for row, key, label, hint in RHS_FIELDS:
            self._add_field(row, key, label, hint)

        # ---- working directory ----
        wd_lbl = tk.Label(self, text="Working Directory = ", anchor="e")
        wd_lbl.grid(row=0, column=3, padx=PADX, pady=PADY, sticky="e")
        self.entry_wd = tk.Entry(self, width=40)
        self.entry_wd.grid(row=0, column=4, padx=PADX, pady=PADY, columnspan=2)
        btn_wd = tk.Button(self, text="Browse", command=self._browse_wd, width=10)
        btn_wd.grid(row=0, column=6, padx=PADX, pady=PADY, sticky="w")

        # ---- output image ----
        self._load_default_image()
        self.label_image = tk.Label(self, text="", image=self._img_result, compound="left")
        self.label_image.image = self._img_result
        self.label_image.grid(row=1, column=3, padx=PADX, pady=PADY, rowspan=16, columnspan=4)

        # ---- status text ----
        self.status_text = tk.Label(self, text="Output Text.........................",
                                    anchor="w")
        self.status_text.grid(row=17, column=3, padx=PADX, pady=PADY, sticky="w", columnspan=4)

        # ---- buttons ----
        btn_load = tk.Button(self, text="LOAD", command=self._load_inputs)
        btn_load.grid(row=18, column=3, padx=PADX, pady=PADY, sticky="w")

        btn_run = tk.Button(self, text="RUN", command=self._run)
        btn_run.grid(row=18, column=4, padx=PADX, pady=PADY, sticky="w")

        self._toggle_var = tk.StringVar(value="Result1")
        toggle_btn = tk.Checkbutton(
            self, text="Toggle images", command=self._toggle_image,
            variable=self._toggle_var, onvalue="Result1", offvalue="Result2",
        )
        toggle_btn.grid(row=18, column=5, padx=PADX, pady=PADY, sticky="w")

        btn_exit = tk.Button(self, text="EXIT", command=self.destroy, width=10)
        btn_exit.grid(row=18, column=6, padx=PADX, pady=PADY, sticky="w")

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

    def _entries_to_params(self) -> GearParams:
        kwargs = {}
        for key, ent in self.entries.items():
            raw = ent.get()
            kwargs[key] = int(raw) if key in _INT_FIELDS else float(raw)
        return GearParams(**kwargs)

    def _read_wd(self) -> str:
        return self.entry_wd.get().strip()

    # ------------------------------------------------------------------ CSV I/O

    def _save_inputs_csv(self) -> None:
        """Write the current UI parameters to work_dir/Inputs.csv."""
        p = self._entries_to_params()
        rows = []
        for key in DEFAULT_PARAMS.__dataclass_fields__:
            rows.append((key, getattr(p, key)))
        df = pd.DataFrame(rows, columns=["parameter", "value"])
        path = os.path.join(self.work_dir, "Inputs.csv")
        df.to_csv(path, index=False)

    def _load_inputs_csv(self) -> GearParams | None:
        """Read work_dir/Inputs.csv and return GearParams, or None."""
        path = os.path.join(self.work_dir, "Inputs.csv")
        if not os.path.exists(path):
            self.status_text.config(text=f"The File is not exists : {path}")
            return None
        df = pd.read_csv(path, index_col="parameter")
        kwargs = {}
        for key in DEFAULT_PARAMS.__dataclass_fields__:
            raw = df.loc[key, "value"]
            kwargs[key] = int(raw) if key in _INT_FIELDS else float(raw)
        self.status_text.config(text=f"The File was loaded : {path}")
        return GearParams(**kwargs)

    # ------------------------------------------------------------------ callbacks

    def _browse_wd(self) -> None:
        wd = filedialog.askdirectory()
        if not wd:
            return
        self.entry_wd.delete(0, "end")
        self.entry_wd.insert(0, wd)
        self.work_dir = wd
        self.status_text.config(text=f"Working Directory : {wd}")

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
        generate(self.params, self.work_dir)
        self._save_inputs_csv()
        self._display_image("Result1.png")
        self.status_text.config(text="Finished generating")

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
        self.status_text.config(text=f"Load Image : {path}")