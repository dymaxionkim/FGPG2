import customtkinter
import os
import numpy as np
import matplotlib.pyplot as plt
import ezdxf
from PIL import Image

PADX = 5
PADY = 2

def button_wd_callback():
    print("button_wd pressed")

def button_load_callback():
    print("button_load pressed")

def button_run_callback():
    print("button_run pressed")

def button_exit_callback():
    print("button_exit pressed")
    exit()

customtkinter.set_default_color_theme("green")
app = customtkinter.CTk()
app.title("FGPG2 with customtkinter")
#app.geometry("640x480")
#app.iconbitmap("./FGPG2.ico")


# Subject
label_x0_1 = customtkinter.CTkLabel(app, text="# Gear Spec", fg_color="transparent", compound="right")
label_x0_1.grid(row=0, column=0, padx=PADX, pady=PADY, sticky="w")

# Module, m
label_m1 = customtkinter.CTkLabel(app, text="Module, m = ", fg_color="transparent", compound="right")
label_m1.grid(row=1, column=0, padx=PADX, pady=PADY, sticky="e")

entry_m = customtkinter.CTkEntry(app, placeholder_text="1.0")
entry_m.grid(row=1, column=1, padx=PADX, pady=PADY)

label_m2 = customtkinter.CTkLabel(app, text="[mm] > 0", fg_color="transparent", compound="left")
label_m2.grid(row=1, column=2, padx=PADX, pady=PADY, sticky="w")

# Teeth Number, z
label_z1 = customtkinter.CTkLabel(app, text="Teeth Number, z = ", fg_color="transparent", compound="right")
label_z1.grid(row=2, column=0, padx=PADX, pady=PADY, sticky="e")

entry_z = customtkinter.CTkEntry(app, placeholder_text="18")
entry_z.grid(row=2, column=1, padx=PADX, pady=PADY)

label_z2 = customtkinter.CTkLabel(app, text="[ea] Minus(-) for Internal Gear", fg_color="transparent", compound="left")
label_z2.grid(row=2, column=2, padx=PADX, pady=PADY, sticky="w")

# Pressure Angle, alpha
label_alpha1 = customtkinter.CTkLabel(app, text="Pressure Angle, alpha = ", fg_color="transparent", compound="right")
label_alpha1.grid(row=3, column=0, padx=PADX, pady=PADY, sticky="e")

entry_alpha = customtkinter.CTkEntry(app, placeholder_text="18")
entry_alpha.grid(row=3, column=1, padx=PADX, pady=PADY)

label_alpha2 = customtkinter.CTkLabel(app, text="[deg] 20 for Standard", fg_color="transparent", compound="left")
label_alpha2.grid(row=3, column=2, padx=PADX, pady=PADY, sticky="w")

# Offset Factor, x
label_x1 = customtkinter.CTkLabel(app, text="Offset Factor, x = ", fg_color="transparent", compound="right")
label_x1.grid(row=4, column=0, padx=PADX, pady=PADY, sticky="e")

entry_x = customtkinter.CTkEntry(app, placeholder_text="0.0")
entry_x.grid(row=4, column=1, padx=PADX, pady=PADY)

label_x2 = customtkinter.CTkLabel(app, text="-1.0 ~ +1.0", fg_color="transparent", compound="left")
label_x2.grid(row=4, column=2, padx=PADX, pady=PADY, sticky="w")

# Backlash Factor, b
label_b1 = customtkinter.CTkLabel(app, text="Backlash Factor, b = ", fg_color="transparent", compound="right")
label_b1.grid(row=5, column=0, padx=PADX, pady=PADY, sticky="e")

entry_b = customtkinter.CTkEntry(app, placeholder_text="0.0")
entry_b.grid(row=5, column=1, padx=PADX, pady=PADY)

label_b2 = customtkinter.CTkLabel(app, text="0.0 ~ +1.0", fg_color="transparent", compound="left")
label_b2.grid(row=5, column=2, padx=PADX, pady=PADY, sticky="w")

# Addendum Factor, a
label_a1 = customtkinter.CTkLabel(app, text="Addendum Factor, a = ", fg_color="transparent", compound="right")
label_a1.grid(row=6, column=0, padx=PADX, pady=PADY, sticky="e")

entry_a = customtkinter.CTkEntry(app, placeholder_text="1.0")
entry_a.grid(row=6, column=1, padx=PADX, pady=PADY)

label_a2 = customtkinter.CTkLabel(app, text="1.0 for Standard", fg_color="transparent", compound="left")
label_a2.grid(row=6, column=2, padx=PADX, pady=PADY, sticky="w")

# Dedendum Factor, d
label_d1 = customtkinter.CTkLabel(app, text="Dedendum Factor, d = ", fg_color="transparent", compound="right")
label_d1.grid(row=7, column=0, padx=PADX, pady=PADY, sticky="e")

entry_d = customtkinter.CTkEntry(app, placeholder_text="1.25")
entry_d.grid(row=7, column=1, padx=PADX, pady=PADY)

label_d2 = customtkinter.CTkLabel(app, text="1.25 for Standard", fg_color="transparent", compound="left")
label_d2.grid(row=7, column=2, padx=PADX, pady=PADY, sticky="w")

# Radius of Hob end, c
label_d1 = customtkinter.CTkLabel(app, text="Radius of Hob end, c = ", fg_color="transparent", compound="right")
label_d1.grid(row=8, column=0, padx=PADX, pady=PADY, sticky="e")

entry_d = customtkinter.CTkEntry(app, placeholder_text="0.2")
entry_d.grid(row=8, column=1, padx=PADX, pady=PADY)

label_d2 = customtkinter.CTkLabel(app, text="[mm]", fg_color="transparent", compound="left")
label_d2.grid(row=8, column=2, padx=PADX, pady=PADY, sticky="w")

# Radius of Tooth end, e
label_d1 = customtkinter.CTkLabel(app, text="Radius of Tooth end, e = ", fg_color="transparent", compound="right")
label_d1.grid(row=9, column=0, padx=PADX, pady=PADY, sticky="e")

entry_d = customtkinter.CTkEntry(app, placeholder_text="0.1")
entry_d.grid(row=9, column=1, padx=PADX, pady=PADY)

label_d2 = customtkinter.CTkLabel(app, text="[mm]", fg_color="transparent", compound="left")
label_d2.grid(row=9, column=2, padx=PADX, pady=PADY, sticky="w")


# Subject
label_x0_1 = customtkinter.CTkLabel(app, text="# Graphics", fg_color="transparent", compound="right")
label_x0_1.grid(row=10, column=0, padx=PADX, pady=PADY, sticky="w")

# Center of Gear , x0
label_x0_1 = customtkinter.CTkLabel(app, text="x0 = ", fg_color="transparent", compound="right")
label_x0_1.grid(row=11, column=0, padx=PADX, pady=PADY, sticky="e")

label_x0 = customtkinter.CTkEntry(app, placeholder_text="0.0")
label_x0.grid(row=11, column=1, padx=PADX, pady=PADY)

label_x0_2 = customtkinter.CTkLabel(app, text="[mm]", fg_color="transparent", compound="left")
label_x0_2.grid(row=11, column=2, padx=PADX, pady=PADY, sticky="w")

# Center of Gear , y0
label_y0_1 = customtkinter.CTkLabel(app, text="y0 = ", fg_color="transparent", compound="right")
label_y0_1.grid(row=12, column=0, padx=PADX, pady=PADY, sticky="e")

entry_y0 = customtkinter.CTkEntry(app, placeholder_text="0.0")
entry_y0.grid(row=12, column=1, padx=PADX, pady=PADY)

label_y0_2 = customtkinter.CTkLabel(app, text="[mm]", fg_color="transparent", compound="left")
label_y0_2.grid(row=12, column=2, padx=PADX, pady=PADY, sticky="w")

# Segmentations , seg_circle
label_seg_circle_1 = customtkinter.CTkLabel(app, text="seg_circle = ", fg_color="transparent", compound="right")
label_seg_circle_1.grid(row=13, column=0, padx=PADX, pady=PADY, sticky="e")

entry_seg_circle = customtkinter.CTkEntry(app, placeholder_text="360")
entry_seg_circle.grid(row=13, column=1, padx=PADX, pady=PADY)

label_seg_circle_2 = customtkinter.CTkLabel(app, text="[ea]", fg_color="transparent", compound="left")
label_seg_circle_2.grid(row=13, column=2, padx=PADX, pady=PADY, sticky="w")

# Segmentations , seg_involute
label_seg_involute_1 = customtkinter.CTkLabel(app, text="seg_involute = ", fg_color="transparent", compound="right")
label_seg_involute_1.grid(row=14, column=0, padx=PADX, pady=PADY, sticky="e")

entry_seg_involute = customtkinter.CTkEntry(app, placeholder_text="15")
entry_seg_involute.grid(row=14, column=1, padx=PADX, pady=PADY)

label_seg_involute_2 = customtkinter.CTkLabel(app, text="[ea]", fg_color="transparent", compound="left")
label_seg_involute_2.grid(row=14, column=2, padx=PADX, pady=PADY, sticky="w")

# Segmentations , seg_edge_r
label_seg_edge_r_1 = customtkinter.CTkLabel(app, text="seg_edge_r = ", fg_color="transparent", compound="right")
label_seg_edge_r_1.grid(row=15, column=0, padx=PADX, pady=PADY, sticky="e")

entry_seg_edge_r = customtkinter.CTkEntry(app, placeholder_text="9")
entry_seg_edge_r.grid(row=15, column=1, padx=PADX, pady=PADY)

label_seg_edge_r_2 = customtkinter.CTkLabel(app, text="[ea]", fg_color="transparent", compound="left")
label_seg_edge_r_2.grid(row=15, column=2, padx=PADX, pady=PADY, sticky="w")

# Segmentations , seg_root_r
label_seg_root_r_1 = customtkinter.CTkLabel(app, text="seg_root_r = ", fg_color="transparent", compound="right")
label_seg_root_r_1.grid(row=16, column=0, padx=PADX, pady=PADY, sticky="e")

entry_seg_root_r = customtkinter.CTkEntry(app, placeholder_text="9")
entry_seg_root_r.grid(row=16, column=1, padx=PADX, pady=PADY)

label_seg_root_r_2 = customtkinter.CTkLabel(app, text="[ea]", fg_color="transparent", compound="left")
label_seg_root_r_2.grid(row=16, column=2, padx=PADX, pady=PADY, sticky="w")

# Segmentations , seg_outer
label_seg_outer_1 = customtkinter.CTkLabel(app, text="seg_outer = ", fg_color="transparent", compound="right")
label_seg_outer_1.grid(row=17, column=0, padx=PADX, pady=PADY, sticky="e")

entry_seg_outer = customtkinter.CTkEntry(app, placeholder_text="5")
entry_seg_outer.grid(row=17, column=1, padx=PADX, pady=PADY)

label_seg_outer_2 = customtkinter.CTkLabel(app, text="[ea]", fg_color="transparent", compound="left")
label_seg_outer_2.grid(row=17, column=2, padx=PADX, pady=PADY, sticky="w")

# Segmentations , seg_root
label_seg_root_1 = customtkinter.CTkLabel(app, text="seg_root = ", fg_color="transparent", compound="right")
label_seg_root_1.grid(row=18, column=0, padx=PADX, pady=PADY, sticky="e")

entry_seg_root = customtkinter.CTkEntry(app, placeholder_text="5")
entry_seg_root.grid(row=18, column=1, padx=PADX, pady=PADY)

label_seg_root_2 = customtkinter.CTkLabel(app, text="[ea]", fg_color="transparent", compound="left")
label_seg_root_2.grid(row=18, column=2, padx=PADX, pady=PADY, sticky="w")

# Scale for one tooth , scale
label_scale_1 = customtkinter.CTkLabel(app, text="seg_root = ", fg_color="transparent", compound="right")
label_scale_1.grid(row=19, column=0, padx=PADX, pady=PADY, sticky="e")

entry_scale = customtkinter.CTkEntry(app, placeholder_text="0.7")
entry_scale.grid(row=19, column=1, padx=PADX, pady=PADY)

label_scale_2 = customtkinter.CTkLabel(app, text="0.1 ~ 1.0", fg_color="transparent", compound="left")
label_scale_2.grid(row=19, column=2, padx=PADX, pady=PADY, sticky="w")


# Subject
label_x0_1 = customtkinter.CTkLabel(app, text="# Output", fg_color="transparent", compound="right")
label_x0_1.grid(row=0, column=3, padx=PADX, pady=PADY, sticky="w")

# Working Directory
label_wd = customtkinter.CTkLabel(app, text="Working Directory = ", fg_color="transparent", compound="right")
label_wd.grid(row=1, column=3, padx=PADX, pady=PADY, sticky="e")

entry_wd = customtkinter.CTkEntry(app, placeholder_text="./Result/")
entry_wd.grid(row=1, column=4, padx=PADX, pady=PADY)

button_wd = customtkinter.CTkButton(app, text="Browse", command=button_wd_callback)
button_wd.grid(row=1, column=5, padx=PADX, pady=PADY, sticky="w")

# Output Image
image_result = customtkinter.CTkImage(light_image=Image.open("./FGPG2.png"), size=(500,500))
label_image = customtkinter.CTkLabel(app, text="", image=image_result, compound="left")
label_image.grid(row=2, column=3, padx=PADX, pady=PADY, rowspan=15, columnspan=3)

# Output Text
label_text = customtkinter.CTkLabel(app, text="Output Text.............................", fg_color="transparent", compound="right")
label_text.grid(row=18, column=3, padx=PADX, pady=PADY, sticky="w", columnspan=3)

# Buttons
button_load = customtkinter.CTkButton(app, text="LOAD", command=button_load_callback)
button_load.grid(row=19, column=3, padx=PADX, pady=PADY, sticky="e")

button_run = customtkinter.CTkButton(app, text="RUN", command=button_run_callback)
button_run.grid(row=19, column=4, padx=PADX, pady=PADY, sticky="ew")

button_exit = customtkinter.CTkButton(app, text="EXIT", command=button_exit_callback)
button_exit.grid(row=19, column=5, padx=PADX, pady=PADY, sticky="w")


app.mainloop()