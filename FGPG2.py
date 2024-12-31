import customtkinter
from customtkinter import filedialog
import os
import numpy as np
import matplotlib.pyplot as plt
import ezdxf
from PIL import Image
import pandas

##############################
"""
# Input Parameters
M,Z,ALPHA,X,B,A,D,C,E,X_0,Y_0,SEG_CIRCLE,SEG_INVOLUTE,SEG_EDGE_R,SEG_ROOT_R,SEG_OUTER,SEG_ROOT,SCALE

# Caculated Parameters
ALPHA_0, ALPHA [rad]
ALPHA_M, Center Line's Slope [Rad]
ALPHA_IS, Start Angle for Involute Curve
THETA_IS, Minimum Range of Parameter to Draw Involute Curve
THETA_IE, Maximum Range of Parameter to Draw Involute Curve
ALPHA_E, Angle between Tooth's Center & End Point of Tooth
X_E, Location of Tooth's End Point
Y_E, Location of Tooth's End Point
X_E0, Location of Edge Round Center
Y_E0, Location of Edge Round Center
THETA3_MIN, Start Angle of Edge Round Curve
THETA3_MAX, End Angle of Edge Round Curve
ALPHA_TS, Start Angle of Root Round Curve
THETA_TE, End Angle of Root Round Curve
P_ANGLE, Gear Pitch Angle
ALIGN_ANGLE, Gear Align Angle

# Linspaces
THETA1, Linspace for Involute curve
THETA3, Linspace for Edge Round curve
THETA_T, Linspace for Root Round Curve
THETA_S, Substitution Variable to plot Root Round Curve
THETA6, Linspace for outer Arc
THETA7, Linspace for root Arc
"""

##############################
# Function
def Internal(Z,X,B,A,D,C,E):
    if Z<0 :
        Z = -Z
        X = -X
        B = -B
        A_temp = A
        D_temp = D
        A = D_temp
        D = A_temp
        C_temp = C
        E_temp = E
        C = E_temp
        E = C_temp
    return Z,X,B,A,D,C,E

def Parameters(M,Z,ALPHA,X,B,A,D,C,E,X_0,Y_0,SEG_CIRCLE,SEG_INVOLUTE,SEG_EDGE_R,SEG_ROOT_R,SEG_OUTER,SEG_ROOT,SCALE):
    ALPHA_0 = ALPHA*(2*np.pi/360)
    ALPHA_M = np.pi/Z
    ALPHA_IS = ALPHA_0+np.pi/(2*Z)+B/(Z*np.cos(ALPHA_0))-(1+2*X/Z)*np.sin(ALPHA_0)/np.cos(ALPHA_0)
    THETA_IS = np.sin(ALPHA_0)/np.cos(ALPHA_0) + 2*(C*(1-np.sin(ALPHA_0))+X-D)/(Z*np.cos(ALPHA_0)*np.sin(ALPHA_0))
    THETA_IE = 2*E/(Z*np.cos(ALPHA_0))+np.sqrt(((Z+2*(X+A-E))/(Z*np.cos(ALPHA_0)))**2-1)
    ALPHA_E = ALPHA_IS + THETA_IE - np.arctan(np.sqrt(((Z+2*(X+A-E))/(Z*np.cos(ALPHA_0)))**2-1))
    X_E = M*((Z/2)+X+A)*np.cos(ALPHA_E)
    Y_E = M*((Z/2)+X+A)*np.sin(ALPHA_E)
    X_E0 = M*(Z/2+X+A-E)*np.cos(ALPHA_E)
    Y_E0 = M*(Z/2+X+A-E)*np.sin(ALPHA_E)
    ALPHA_TS = (2*(C*(1-np.sin(ALPHA_0))-D)*np.sin(ALPHA_0)+B)/(Z*np.cos(ALPHA_0))-2*C*np.cos(ALPHA_0)/Z+np.pi/(2*Z)
    THETA_TE = 2*C*np.cos(ALPHA_0)/Z - 2*(D-X-C*(1-np.sin(ALPHA_0)))*np.cos(ALPHA_0)/(Z*np.sin(ALPHA_0))
    # modify "E"
    if (ALPHA_E>ALPHA_M) and (ALPHA_M>ALPHA_IS+THETA_IE-np.arctan(THETA_IE)) :
        E = (E/2)*np.cos(ALPHA_0)*(THETA_IE-np.sqrt((1/np.cos(ALPHA_IS+THETA_IE-ALPHA_M))**2-1))
    P_ANGLE = 2*np.pi/Z
    ALIGN_ANGLE = np.pi/2-np.pi/Z
    return ALPHA_0,ALPHA_M,ALPHA_IS,THETA_IS,THETA_IE,ALPHA_E,X_E,Y_E,X_E0,Y_E0,ALPHA_TS,THETA_TE,E,P_ANGLE,ALIGN_ANGLE

def SymetryY(XX,YY):
    XX2 =  XX[::-1]
    YY2 = -YY[::-1]
    return XX2,YY2

def InvoluteCurve(M,Z,SEG_INVOLUTE,THETA_IS,THETA_IE,ALPHA_0,ALPHA_IS):
    THETA1 = np.linspace(THETA_IS,THETA_IE,SEG_INVOLUTE)
    X11 = np.ones(len(THETA1))
    X11 = (1/2)*M*Z*np.cos(ALPHA_0)*np.sqrt(1+THETA1**2)*np.cos(ALPHA_IS+THETA1-np.arctan(THETA1))
    Y11 = (1/2)*M*Z*np.cos(ALPHA_0)*np.sqrt(1+THETA1**2)*np.sin(ALPHA_IS+THETA1-np.arctan(THETA1))
    return X11,Y11,THETA1

def EdgeRoundCurve(M,E,X11,Y11,X_E,Y_E,X_E0,Y_E0,SEG_EDGE_R):
    THETA3_MIN = np.arctan((Y11[len(Y11)-1]-Y_E0)/(X11[len(X11)-1]-X_E0))
    THETA3_MAX = np.arctan((Y_E-Y_E0)/(X_E-X_E0))
    THETA3 = np.linspace(THETA3_MIN,THETA3_MAX,SEG_EDGE_R)
    X21 = M*E*np.cos(THETA3) + X_E0
    Y21 = M*E*np.sin(THETA3) + Y_E0
    return X21,Y21,THETA3,THETA3_MIN,THETA3_MAX

def RootRoundCurve(M,Z,X,D,C,B,THETA_TE,ALPHA_TS,SEG_ROOT_R):
    THETA_T = np.linspace(0,THETA_TE,SEG_ROOT_R)
    if (C!=0) and ((D-X-C)==0) :
        # mc瑜� 諛섏��由꾩쑝濡� �븯�뒗 �썝�샇瑜� 洹몃젮�꽌 ���泥댄븯寃� �맖
        THETA_S = (np.pi/2)*np.ones(len(THETA_T))
    elif (D-X-C)!=0 :
        THETA_S = np.arctan((M*Z*THETA_T/2)/(M*D-M*X-M*C))
    X31 = M*((Z/2+X-D+C)*np.cos(THETA_T+ALPHA_TS)+(Z/2)*THETA_T*np.sin(THETA_T+ALPHA_TS)-C*np.cos(THETA_S+THETA_T+ALPHA_TS))
    Y31 = M*((Z/2+X-D+C)*np.sin(THETA_T+ALPHA_TS)-(Z/2)*THETA_T*np.cos(THETA_T+ALPHA_TS)-C*np.sin(THETA_S+THETA_T+ALPHA_TS))
    return X31,Y31,THETA_T,THETA_S

def OuterArc(M,Z,X,A,ALPHA_E,ALPHA_M,SEG_OUTER):
    THETA6 = np.linspace(ALPHA_E,ALPHA_M,SEG_OUTER) 
    X41 = M*(Z/2+A+X)*np.cos(THETA6)
    Y41 = M*(Z/2+A+X)*np.sin(THETA6)
    return X41,Y41,THETA6

def RootArc(M,Z,X,D,ALPHA_TS,SEG_ROOT):
    THETA7 = np.linspace(0,ALPHA_TS,SEG_ROOT) 
    X51 = M*(Z/2-D+X)*np.cos(THETA7)
    Y51 = M*(Z/2-D+X)*np.sin(THETA7)
    return X51,Y51,THETA7

def CombineTooth(X11,Y11,X21,Y21,X31,Y31,X41,Y41,X51,Y51,X12,Y12,X22,Y22,X32,Y32,X42,Y42,X52,Y52):
    # Combine to one tooth
    X23 = np.delete(X22,0)
    Y23 = np.delete(Y22,0)
    X13 = np.delete(X12,0)
    Y13 = np.delete(Y12,0)
    X33 = np.delete(X32,0)
    Y33 = np.delete(Y32,0)
    X34 = np.delete(X31,0)
    Y34 = np.delete(Y31,0)
    X14 = np.delete(X11,0)
    Y14 = np.delete(Y11,0)
    X24 = np.delete(X21,0)
    Y24 = np.delete(Y21,0)
    X1 = np.concatenate((X42,X23,X13,X33,X52,X51,X34,X14,X24,X41))
    Y1 = np.concatenate((Y42,Y23,Y13,Y33,Y52,Y51,Y34,Y14,Y24,Y41))
    return X1,Y1

def Transform(Xtemp,Ytemp,X_0,Y_0):
    Xtemp = Xtemp + X_0
    Ytemp = Ytemp + Y_0
    return Xtemp,Ytemp

def Rotation(Xtemp,Ytemp,ANGLE,i):
    XX = np.cos(ANGLE*i)*Xtemp - np.sin(ANGLE*i)*Ytemp
    YY = np.sin(ANGLE*i)*Xtemp + np.cos(ANGLE*i)*Ytemp
    return XX,YY

def Circle(DIA,SEG_CIRCLE):
    THETA0 = np.linspace(0.0,2*np.pi,SEG_CIRCLE)
    XX = DIA/2*np.sin(THETA0)
    YY = DIA/2*np.cos(THETA0)
    return XX,YY

def CombineToothDXF(X12,Y12,X22,Y22,X32,Y32):
    X13 = np.delete(X12,0)
    Y13 = np.delete(Y12,0)
    X33 = np.delete(X32,0)
    Y33 = np.delete(Y32,0)
    X1 = np.concatenate((X22,X13,X33))
    Y1 = np.concatenate((Y22,Y13,Y33))
    return X1,Y1

def SaveDXF(X12,Y12,X22,Y22,X32,Y32,Z,X_0,Y_0,OUTER_DIA,ROOT_DIA,OFFSET_DIA):
    doc = ezdxf.new('R2000')
    msp = doc.modelspace()
    # Spline Curve
    Xdxf1,Ydxf1 = CombineToothDXF(X12,Y12,X22,Y22,X32,Y32)
    ALIGN_ANGLE_DXF = np.pi/2+np.pi/Z
    Xdxf3,Ydxf3 = Rotation(Xdxf1,Ydxf1,ALIGN_ANGLE_DXF,1) # Align to top
    Xdxf4 = -Xdxf3
    Ydxf4 =  Ydxf3
    Xdxf5,Ydxf5 = Transform(Xdxf3,Ydxf3,X_0,Y_0)
    Xdxf6,Ydxf6 = Transform(Xdxf4,Ydxf4,X_0,Y_0)
    cpoint5 = [(Xdxf5[0],Ydxf5[0])]
    cpoint6 = [(Xdxf6[0],Ydxf6[0])]
    for i in range(0,len(Xdxf5)) :
        cpoint5.append((Xdxf5[i],Ydxf5[i]))
        cpoint6.append((Xdxf6[i],Ydxf6[i]))
    msp.add_spline(cpoint5)
    msp.add_spline(cpoint6)
    # Arc outer
    outer_angle = np.arctan(Xdxf3[0]/Ydxf3[0])*360.0/(2.0*np.pi)
    msp.add_arc(center=(X_0,Y_0), radius=OUTER_DIA/2, start_angle=90+outer_angle, end_angle=90-outer_angle)
    # Arc root
    root_start_angle = 180/Z
    root_end_angle = np.arctan(Xdxf3[len(Xdxf3)-1]/Ydxf3[len(Xdxf3)-1])*360.0/(2.0*np.pi)
    msp.add_arc(center=(X_0,Y_0), radius=ROOT_DIA/2, start_angle=90-root_start_angle, end_angle=90+root_end_angle)
    msp.add_arc(center=(X_0,Y_0), radius=ROOT_DIA/2, start_angle=90-root_end_angle, end_angle=90+root_start_angle)
    # Line dividing
    Xdivide1 = -np.sin(np.pi/Z)*ROOT_DIA/2 + X_0
    Ydivide1 =  np.cos(np.pi/Z)*ROOT_DIA/2 + Y_0
    Xdivide2 = -np.sin(-np.pi/Z)*ROOT_DIA/2 + X_0
    Ydivide2 =  np.cos(-np.pi/Z)*ROOT_DIA/2 + Y_0
    msp.add_line([X_0,Y_0],[Xdivide1,Ydivide1])
    msp.add_line([X_0,Y_0],[Xdivide2,Ydivide2])
    # Circle offset
    msp.add_circle((X_0,Y_0),radius=OFFSET_DIA/2)
    # Output
    Result = os.path.join(WorkingDirectory, f'Result.dxf')
    doc.saveas(Result)

def SaveSpec(M,ALPHA,Z,X,A,D,B,BASE_DIA,PITCH_DIA,OFFSET_DIA,ROOT_DIA,OUTER_DIA):
    Result3 = os.path.join(WorkingDirectory, f'Result.csv')
    fileout = open(Result3, "w")
    if ALPHA == 20 :
        fileout.write("Type,"+"Standard,\n")
    else :
        fileout.write("Type,"+"Non-Standard,\n")
    fileout.write("Module,"+repr(M)+",mm"+"\n")
    fileout.write("Pressure Angle,"+repr(ALPHA)+",deg\n")
    fileout.write("Teeth Number,"+repr(Z)+",ea\n")
    fileout.write("Offset Factor,"+repr(X)+",\n")
    fileout.write("Offset,"+repr(X*M)+",mm\n")
    fileout.write("Backlash Factor,"+repr(B)+",\n")
    fileout.write("Backlash,"+repr(B*M)+",mm\n")
    fileout.write("Addendum Factor,"+repr(A)+",\n")
    fileout.write("Addendum,"+repr(A*M)+",mm\n")
    fileout.write("Dedendum Factor,"+repr(D)+",\n")
    fileout.write("Dedendum,"+repr(D*M)+",mm\n")
    fileout.write("Total Tooth Height,"+repr(A*M+D*M)+",mm\n")
    fileout.write("Base Circle Dia,"+str(BASE_DIA)+",mm\n")
    fileout.write("Pitch Circle Dia,"+repr(PITCH_DIA)+",mm\n")
    fileout.write("Offset Circle Dia,"+repr(OFFSET_DIA)+",mm\n")
    fileout.write("Root Circle Dia,"+repr(ROOT_DIA)+",mm\n")
    fileout.write("Outer Circle Dia,"+repr(OUTER_DIA)+",mm\n")
    fileout.close()

def FGPG2_PLOT(M,Z,ALPHA,X,B,A,D,C,E,X_0,Y_0,SEG_CIRCLE,SEG_INVOLUTE,SEG_EDGE_R,SEG_ROOT_R,SEG_OUTER,SEG_ROOT,SCALE):

    ####################
    # Set Graphics
    fig = plt.figure(figsize=(5,5))
    plt.axes().set_aspect('equal')
    plt.title('Fine Gear Profile Generator 2')
    plt.grid(True)
    
    # Gear tooth
    Z,X,B,A,D,C,E = Internal(Z,X,B,A,D,C,E)
    ALPHA_0,ALPHA_M,ALPHA_IS,THETA_IS,THETA_IE,ALPHA_E,X_E,Y_E,X_E0,Y_E0,ALPHA_TS,THETA_TE,E,P_ANGLE,ALIGN_ANGLE = Parameters(M,Z,ALPHA,X,B,A,D,C,E,X_0,Y_0,SEG_CIRCLE,SEG_INVOLUTE,SEG_EDGE_R,SEG_ROOT_R,SEG_OUTER,SEG_ROOT,SCALE)
    X11,Y11,THETA1 = InvoluteCurve(M,Z,SEG_INVOLUTE,THETA_IS,THETA_IE,ALPHA_0,ALPHA_IS)
    X12,Y12 = SymetryY(X11,Y11)
    X21,Y21,THETA3,THETA3_MIN,THETA3_MAX = EdgeRoundCurve(M,E,X11,Y11,X_E,Y_E,X_E0,Y_E0,SEG_EDGE_R)
    X22,Y22 = SymetryY(X21,Y21)
    X31,Y31,THETA_T,THETA_S = RootRoundCurve(M,Z,X,D,C,B,THETA_TE,ALPHA_TS,SEG_ROOT_R)
    X32,Y32 = SymetryY(X31,Y31)
    X41,Y41,THETA6 = OuterArc(M,Z,X,A,ALPHA_E,ALPHA_M,SEG_OUTER)
    X42,Y42 = SymetryY(X41,Y41)
    X51,Y51,THETA7 = RootArc(M,Z,X,D,ALPHA_TS,SEG_ROOT)
    X52,Y52 = SymetryY(X51,Y51)

    # Plot Whole Gear
    X1,Y1 = CombineTooth(X11,Y11,X21,Y21,X31,Y31,X41,Y41,X51,Y51,X12,Y12,X22,Y22,X32,Y32,X42,Y42,X52,Y52)
    ALIGN_ANGLE = np.pi/2-np.pi/Z
    X2,Y2 = Rotation(X1,Y1,ALIGN_ANGLE,1) # Align to top
    # Array
    P_ANGLE = 2*np.pi/Z
    for i in range(0,Z):
        Xtemp,Ytemp = Rotation(X2,Y2,P_ANGLE,i)
        Xtemp,Ytemp = Transform(Xtemp,Ytemp,X_0,Y_0)
        plt.plot(Xtemp,Ytemp,'-',linewidth=1.5,color='black',label='_nolegend_')
    
    # Base Cicle
    base_dia = M*Z*np.cos(ALPHA_0)
    X_base,Y_base = Circle(base_dia,SEG_CIRCLE)
    X_base,Y_base = Transform(X_base,Y_base,X_0,Y_0)
    plt.plot(X_base, Y_base, ':', linewidth=1.0, color='cyan',label='Base Cicle')
    
    # Pitch Circle
    pitch_dia = M*Z
    X_pitch,Y_pitch = Circle(pitch_dia,SEG_CIRCLE)
    X_pitch,Y_pitch = Transform(X_pitch,Y_pitch,X_0,Y_0)
    plt.plot(X_pitch, Y_pitch, ':', linewidth=1.0, color='magenta',label='Pitch Circle')

    # Offset Circle
    offset_dia = 2*M*(Z/2+X)
    X_offset,Y_offset = Circle(offset_dia,SEG_CIRCLE)
    X_offset,Y_offset = Transform(X_offset,Y_offset,X_0,Y_0)
    plt.plot(X_offset, Y_offset, ':', linewidth=1.0, color='red',label='Offset Circle')

    # Outer Circle
    outer_dia = 2*M*(Z/2+X+A)
    X_out,Y_out = Circle(outer_dia,SEG_CIRCLE)
    X_out,Y_out = Transform(X_out,Y_out,X_0,Y_0)
    plt.plot(X_out, Y_out, ':', linewidth=1.0, color='brown',label='Outer Circle')

    # Root Circle
    root_dia = 2*M*(Z/2+X-D)
    X_root,Y_root = Circle(root_dia,SEG_CIRCLE)
    X_root,Y_root = Transform(X_root,Y_root,X_0,Y_0)
    plt.plot(X_root, Y_root, ':', linewidth=1.0, color='grey',label='Root Circle')
    
    # Annotate
    Cheight = root_dia/22.5
    Nrow = 15*Cheight
    plt.text(X_0,Y_0+Nrow/2-Cheight*0, 'Module m=%s[mm]'%(M),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*1, 'Teeth Number z=%s[ea]'%(Z),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*2, 'Pressure angle alpha=%s[deg]'%(ALPHA),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*3, 'Offset factor x=%s'%(X),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*4, 'Backlash factor b=%s'%(B),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*5, 'Addendum factor a=%s'%(A),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*6, 'Dedendum factor d=%s'%(D),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*7, 'Radius Factor of Edge Round of Hob c=%s'%(C),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*8, 'Radius Factor of Edge Round of Tooth e=%s'%(E),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*9, 'Center Position = %s[mm]'%([X_0,Y_0]),
        verticalalignment='center', horizontalalignment='center', color='green', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*10, 'Base Circle Dia = %s[mm]'%(base_dia),
        verticalalignment='center', horizontalalignment='center', color='cyan', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*11, 'Pitch Circle Dia = %s[mm]'%(pitch_dia),
        verticalalignment='center', horizontalalignment='center', color='magenta', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*12, 'Offset Circle Dia = %s[mm]'%(offset_dia),
        verticalalignment='center', horizontalalignment='center', color='red', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*13, 'Outer Circle Dia = %s[mm]'%(outer_dia),
        verticalalignment='center', horizontalalignment='center', color='brown', fontsize="x-small")
    plt.text(X_0,Y_0+Nrow/2-Cheight*14, 'Root Circle Dia = %s[mm]'%(root_dia),
        verticalalignment='center', horizontalalignment='center', color='grey', fontsize="x-small")
            
    # Save Figure for Total Gear
    Result = os.path.join(WorkingDirectory, f'Result1.png')
    plt.savefig(Result,dpi=100)
    #Result = os.path.join(WorkingDirectory, f'Result1.svg')
    #plt.savefig(Result)

    # Save Figure for One Tooth
    Y_height = ((outer_dia-root_dia)/2.0)/SCALE
    Y_cen = (outer_dia+root_dia)/4.0 + Y_0
    Y_min = Y_cen - Y_height/2.0
    Y_max = Y_cen + Y_height/2.0
    X_min = -(Y_height/2.0)+X_0
    X_max = (Y_height/2.0)+X_0
    plt.xlim(X_min, X_max)
    plt.ylim(Y_min, Y_max)
    Result2 = os.path.join(WorkingDirectory, f'Result2.png')
    plt.savefig(Result2,dpi=100)
    #Result2 = os.path.join(WorkingDirectory, f'Result2.svg')
    #plt.savefig(Result2)

    ####################
    # Save DXF
    SaveDXF(X12,Y12,X22,Y22,X32,Y32,Z,X_0,Y_0,outer_dia,root_dia,offset_dia)
    
    ####################
    # Save Gear Spec
    SaveSpec(M,ALPHA,Z,X,A,D,B,base_dia,pitch_dia,offset_dia,root_dia,outer_dia)


##############################

# Gap between pads in customtkinter
PADX = 1
PADY = 1

# Functions
def init_parameters():
    entry_m.insert(0,1.0)
    entry_z.insert(0,18)
    entry_alpha.insert(0,20)
    entry_x.insert(0,0.0)
    entry_b.insert(0,0.05)
    entry_a.insert(0,1.0)
    entry_d.insert(0,1.25)
    entry_c.insert(0,0.2)
    entry_e.insert(0,0.1)
    entry_x0.insert(0,0.0)
    entry_y0.insert(0,0.0)
    entry_seg_circle.insert(0,360)
    entry_seg_involute.insert(0,15)
    entry_seg_edge_r.insert(0,5)
    entry_seg_root_r.insert(0,5)
    entry_seg_outer.insert(0,5)
    entry_seg_root.insert(0,5)
    entry_scale.insert(0,0.7)
    entry_wd.insert(0,"./Result")

def read_parameters():
    global m, z, alpha, x, b, a, d, c, e, x_0, y_0, seg_circle, seg_involute, seg_edge_r, seg_root_r, seg_outer, seg_root, scale, WorkingDirectory
    m = float(entry_m.get())
    z = int(entry_z.get())
    alpha = float(entry_alpha.get())
    x = float(entry_x.get())
    b = float(entry_b.get())
    a = float(entry_a.get())
    d = float(entry_d.get())
    c = float(entry_c.get())
    e = float(entry_e.get())
    x_0 = float(entry_x0.get())
    y_0 = float(entry_y0.get())
    seg_circle = int(entry_seg_circle.get())
    seg_involute = int(entry_seg_involute.get())
    seg_edge_r = int(entry_seg_edge_r.get())
    seg_root_r = int(entry_seg_root_r.get())
    seg_outer = int(entry_seg_outer.get())
    seg_root = int(entry_seg_root.get())
    scale = float(entry_scale.get())
    WorkingDirectory = entry_wd.get()

def load_parameters():
    global m, z, alpha, x, b, a, d, c, e, x_0, y_0, seg_circle, seg_involute, seg_edge_r, seg_root_r, seg_outer, seg_root, scale, Inputs
    parameters = pandas.read_csv(Inputs,index_col="parameter")
    m = parameters.loc['m'].astype('float').values[0]
    z = parameters.loc['z'].astype('int').values[0]
    alpha = parameters.loc['alpha'].astype('float').values[0]
    x = parameters.loc['x'].astype('float').values[0]
    b = parameters.loc['b'].astype('float').values[0]
    a = parameters.loc['a'].astype('float').values[0]
    d = parameters.loc['d'].astype('float').values[0]
    c = parameters.loc['c'].astype('float').values[0]
    e = parameters.loc['e'].astype('float').values[0]
    x_0 = parameters.loc['x_0'].astype('float').values[0]
    y_0 = parameters.loc['y_0'].astype('float').values[0]
    seg_circle = parameters.loc['seg_circle'].astype('int').values[0]
    seg_involute = parameters.loc['seg_involute'].astype('int').values[0]
    seg_edge_r = parameters.loc['seg_edge_r'].astype('int').values[0]
    seg_root_r = parameters.loc['seg_root_r'].astype('int').values[0]
    seg_outer = parameters.loc['seg_outer'].astype('int').values[0]
    seg_root = parameters.loc['seg_root'].astype('int').values[0]
    scale = parameters.loc['scale'].astype('float').values[0]
    print(parameters)
    # Apply into UI
    entry_m.delete(0,last_index='end')
    entry_m.insert(0,m)
    entry_z.delete(0,last_index='end')
    entry_z.insert(0,z)
    entry_alpha.delete(0,last_index='end')
    entry_alpha.insert(0,alpha)
    entry_x.delete(0,last_index='end')
    entry_x.insert(0,x)
    entry_b.delete(0,last_index='end')
    entry_b.insert(0,b)
    entry_a.delete(0,last_index='end')
    entry_a.insert(0,a)
    entry_d.delete(0,last_index='end')
    entry_d.insert(0,d)
    entry_c.delete(0,last_index='end')
    entry_c.insert(0,c)
    entry_e.delete(0,last_index='end')
    entry_e.insert(0,e)
    entry_x0.delete(0,last_index='end')
    entry_x0.insert(0,x_0)
    entry_y0.delete(0,last_index='end')
    entry_y0.insert(0,y_0)
    entry_seg_circle.delete(0,last_index='end')
    entry_seg_circle.insert(0,seg_circle)
    entry_seg_involute.delete(0,last_index='end')
    entry_seg_involute.insert(0,seg_involute)
    entry_seg_edge_r.delete(0,last_index='end')
    entry_seg_edge_r.insert(0,seg_edge_r)
    entry_seg_root_r.delete(0,last_index='end')
    entry_seg_root_r.insert(0,seg_root_r)
    entry_seg_outer.delete(0,last_index='end')
    entry_seg_outer.insert(0,seg_outer)
    entry_seg_root.delete(0,last_index='end')
    entry_seg_root.insert(0,seg_root)
    entry_scale.delete(0,last_index='end')
    entry_scale.insert(0,scale)

def save_parameters():
    global m, z, alpha, x, b, a, d, c, e, x_0, y_0, seg_circle, seg_involute, seg_edge_r, seg_root_r, seg_outer, seg_root, scale, WorkingDirectory
    read_parameters()
    Outputs = os.path.join(WorkingDirectory, f'Inputs.csv')
    parameters = pandas.read_csv('parameters.csv',index_col="parameter")
    parameters.loc['m','value'] = m
    parameters.loc['z','value'] = z
    parameters.loc['alpha','value'] = alpha
    parameters.loc['x','value'] = x
    parameters.loc['b','value'] = b
    parameters.loc['a','value'] = a
    parameters.loc['d','value'] = d
    parameters.loc['c','value'] = c
    parameters.loc['e','value'] = e
    parameters.loc['x_0','value'] = x_0
    parameters.loc['y_0','value'] = y_0
    parameters.loc['seg_circle','value'] = seg_circle
    parameters.loc['seg_involute','value'] = seg_involute
    parameters.loc['seg_edge_r','value'] = seg_edge_r
    parameters.loc['seg_root_r','value'] = seg_root_r
    parameters.loc['seg_outer','value'] = seg_outer
    parameters.loc['seg_root','value'] = seg_root
    parameters.loc['scale','value'] = scale
    parameters.to_csv(Outputs)
    print(parameters)

# Callback Functions
def button_wd_callback():
    print("button_wd pressed")
    global WorkingDirectory
    WorkingDirectory = filedialog.askdirectory()
    entry_wd.delete(0,last_index='end')
    entry_wd.insert(0,WorkingDirectory)
    print('Working Directory : %s'%WorkingDirectory)
    label_text.configure(text='Working Directory : %s'%WorkingDirectory)

def button_load_callback():
    print("button_load pressed")
    global Inputs
    read_parameters()
    Inputs = os.path.join(WorkingDirectory, f'Inputs.csv')
    if os.path.exists(Inputs) :
        load_parameters()
        label_text.configure(text='The File was loaded : %s'%Inputs)
    else :
        label_text.configure(text='The File is not exists : %s'%Inputs)

def button_run_callback():
    global label_image
    print("button_run pressed")
    read_parameters()
    os.makedirs(WorkingDirectory, exist_ok=True)
    FGPG2_PLOT(m,z,alpha,x,b,a,d,c,e,x_0,y_0,seg_circle,seg_involute,seg_edge_r,seg_root_r,seg_outer,seg_root,scale)
    save_parameters()
    Result = os.path.join(WorkingDirectory, f'Result1.png')
    image_result = customtkinter.CTkImage(light_image=Image.open(Result), size=(500,500))
    label_image = customtkinter.CTkLabel(app, text="", image=image_result)
    label_image.grid(row=1, column=3, padx=PADX, pady=PADY, rowspan=16, columnspan=4)
    label_text.configure(text='Finished geneating')

def switch_event():
    global WorkingDirectory
    print("switch toggled, current value:", switch_toggle_var.get())
    if switch_toggle_var.get()=="Result1" :
        Result = os.path.join(WorkingDirectory, f'Result1.png')
    else :
        Result = os.path.join(WorkingDirectory, f'Result2.png')
    image_result = customtkinter.CTkImage(light_image=Image.open(Result), size=(500,500))
    label_image = customtkinter.CTkLabel(app, text="", image=image_result)
    label_image.grid(row=1, column=3, padx=PADX, pady=PADY, rowspan=16, columnspan=4)
    label_text.configure(text='load Image : %s'%Result)

def button_exit_callback():
    print("button_exit pressed")
    exit()

################
# GUI
customtkinter.set_default_color_theme("green")
app = customtkinter.CTk()
app.title("FGPG2 with customtkinter")
app.geometry("915x635")
app.resizable(width=False, height=False)
app.iconbitmap('FGPG2.ico')
font16 = customtkinter.CTkFont(size=16)

# Subject
label_x0_1 = customtkinter.CTkLabel(app, text="# Gear Spec", fg_color="transparent", compound="right", font=font16)
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

label_z2 = customtkinter.CTkLabel(app, text="[ea] : - for Internal", fg_color="transparent", compound="left")
label_z2.grid(row=2, column=2, padx=PADX, pady=PADY, sticky="w")

# Pressure Angle, alpha
label_alpha1 = customtkinter.CTkLabel(app, text="Pressure Angle, alpha = ", fg_color="transparent", compound="right")
label_alpha1.grid(row=3, column=0, padx=PADX, pady=PADY, sticky="e")

entry_alpha = customtkinter.CTkEntry(app, placeholder_text="18")
entry_alpha.grid(row=3, column=1, padx=PADX, pady=PADY)

label_alpha2 = customtkinter.CTkLabel(app, text="[deg] : 20 for Standard", fg_color="transparent", compound="left")
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
label_c1 = customtkinter.CTkLabel(app, text="Radius of Hob end, c = ", fg_color="transparent", compound="right")
label_c1.grid(row=8, column=0, padx=PADX, pady=PADY, sticky="e")

entry_c = customtkinter.CTkEntry(app, placeholder_text="0.2")
entry_c.grid(row=8, column=1, padx=PADX, pady=PADY)

label_c2 = customtkinter.CTkLabel(app, text="[mm]", fg_color="transparent", compound="left")
label_c2.grid(row=8, column=2, padx=PADX, pady=PADY, sticky="w")

# Radius of Tooth end, e
label_e1 = customtkinter.CTkLabel(app, text="Radius of Tooth end, e = ", fg_color="transparent", compound="right")
label_e1.grid(row=9, column=0, padx=PADX, pady=PADY, sticky="e")

entry_e = customtkinter.CTkEntry(app, placeholder_text="0.1")
entry_e.grid(row=9, column=1, padx=PADX, pady=PADY)

label_e2 = customtkinter.CTkLabel(app, text="[mm]", fg_color="transparent", compound="left")
label_e2.grid(row=9, column=2, padx=PADX, pady=PADY, sticky="w")


# Subject
label_x0_1 = customtkinter.CTkLabel(app, text="# Graphics", fg_color="transparent", compound="right", font=font16)
label_x0_1.grid(row=10, column=0, padx=PADX, pady=PADY, sticky="w")

# Center of Gear , x0
label_x0_1 = customtkinter.CTkLabel(app, text="x0 = ", fg_color="transparent", compound="right")
label_x0_1.grid(row=11, column=0, padx=PADX, pady=PADY, sticky="e")

entry_x0 = customtkinter.CTkEntry(app, placeholder_text="0.0")
entry_x0.grid(row=11, column=1, padx=PADX, pady=PADY)

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
label_scale_1 = customtkinter.CTkLabel(app, text="scale = ", fg_color="transparent", compound="right")
label_scale_1.grid(row=19, column=0, padx=PADX, pady=PADY, sticky="e")

entry_scale = customtkinter.CTkEntry(app, placeholder_text="0.7")
entry_scale.grid(row=19, column=1, padx=PADX, pady=PADY)

label_scale_2 = customtkinter.CTkLabel(app, text="0.1 ~ 1.0", fg_color="transparent", compound="left")
label_scale_2.grid(row=19, column=2, padx=PADX, pady=PADY, sticky="w")


# Working Directory
label_wd = customtkinter.CTkLabel(app, text="Working Directory = ", fg_color="transparent", compound="right")
label_wd.grid(row=0, column=3, padx=PADX, pady=PADY, sticky="e")

entry_wd = customtkinter.CTkEntry(app, placeholder_text="./Result/", width=300)
entry_wd.grid(row=0, column=4, padx=PADX, pady=PADY, columnspan=2)

button_wd = customtkinter.CTkButton(app, text="Browse", command=button_wd_callback, width=50)
button_wd.grid(row=0, column=6, padx=PADX, pady=PADY, sticky="w")

# Output Image
image_result = customtkinter.CTkImage(light_image=Image.open("./FGPG2.png"), size=(500,500))
label_image = customtkinter.CTkLabel(app, text="", image=image_result, compound="left")
label_image.grid(row=1, column=3, padx=PADX, pady=PADY, rowspan=16, columnspan=4)

# Output Text
label_text = customtkinter.CTkLabel(app, text="Output Text.............................", fg_color="transparent", compound="right")
label_text.grid(row=17, column=3, padx=PADX, pady=PADY, sticky="w", columnspan=4)

# Buttons
button_load = customtkinter.CTkButton(app, text="LOAD", command=button_load_callback)
button_load.grid(row=18, column=3, padx=PADX, pady=PADY, sticky="w")

button_run = customtkinter.CTkButton(app, text="RUN", command=button_run_callback)
button_run.grid(row=18, column=4, padx=PADX, pady=PADY, sticky="w")

switch_toggle_var = customtkinter.StringVar(value="Result1")
switch_toggle = customtkinter.CTkSwitch(app, text="Toggle images", command=switch_event, variable=switch_toggle_var, onvalue="Result1", offvalue="Result2")
switch_toggle.grid(row=18, column=5, padx=PADX, pady=PADY, sticky="w")

button_exit = customtkinter.CTkButton(app, text="EXIT", command=button_exit_callback, width=50)
button_exit.grid(row=18, column=6, padx=PADX, pady=PADY, sticky="w")

# Init
init_parameters()
read_parameters()

app.mainloop()