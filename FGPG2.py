import os
import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
import ezdxf

##############################
# Function
def FGPG2_PLOT(M,Z,ALPHA,X,B,A,D,C,E,X_0,Y_0,SEG_CIRCLE,SEG_INVOLUTE,SEG_EDGE_R,SEG_ROOT_R,SEG_CENTER,SEG_OUTER,SEG_ROOT,SCALE) :

    ####################
    # Graphics
    fig = plt.figure(figsize=(5,5))
    plt.axes().set_aspect('equal')
    plt.title('Fine Gear Profile Generator 2')
    plt.grid('On')

    alpha_0 = ALPHA*(2*np.pi/360)

    # Involute Curve
    # alpha_m = Center Line's Slope [Rad]
    alpha_m = np.pi/Z
    # alpha_is = Start Angle for Involute Curve
    alpha_is = alpha_0 + np.pi/(2*Z) + B/(Z*np.cos(alpha_0)) - (1+2*X/Z)*np.sin(alpha_0)/np.cos(alpha_0)
    # theta_is = Minimum Range of Parameter to Draw Involute Curve
    theta_is = np.sin(alpha_0)/np.cos(alpha_0) + 2*(C*(1-np.sin(alpha_0))+X-D)/(Z*np.cos(alpha_0)*np.sin(alpha_0))
    # theta_ie = Maximum Range of Parameter to Draw Involute Curve
    theta_ie = 2*E/(Z*np.cos(alpha_0)) + np.sqrt( ((Z+2*(X+A-E))/(Z*np.cos(alpha_0)))**2 - 1 )
    THETA1 = np.linspace(theta_is,theta_ie,SEG_INVOLUTE) 
    X11 = np.ones(len(THETA1))
    X11 = (1/2)*M*Z*np.cos(alpha_0)*np.sqrt(1+THETA1**2)*np.cos(alpha_is+THETA1-np.arctan(THETA1))
    Y11 = (1/2)*M*Z*np.cos(alpha_0)*np.sqrt(1+THETA1**2)*np.sin(alpha_is+THETA1-np.arctan(THETA1))
    X12 = X11
    Y12 = -Y11
    X12 = X12[::-1]
    Y12 = Y12[::-1]

    # Edge Round Curve of Tooth
    # alpha_e = Angle between Tooth's Center & End Point of Tooth
    alpha_e = alpha_is + theta_ie - np.arctan(np.sqrt( ( (Z+2*(X+A-E))/(Z*np.cos(alpha_0)) )**2 - 1 ))
    # modify "E"
    if (alpha_e>alpha_m) and (alpha_m>alpha_is+theta_ie-np.arctan(theta_ie)) :
        E = (E/2)*np.cos(alpha_0)*( theta_ie - np.sqrt( (1/np.cos(alpha_is+theta_ie-alpha_m))**2-1 ) )
    # x_e, y_e = Location of Tooth's End Point
    x_e = M*((Z/2)+X+A)*np.cos(alpha_e)
    y_e = M*((Z/2)+X+A)*np.sin(alpha_e)
    # x_e0, y_e0 = Location of Edge Round Center
    x_e0 = M*(Z/2+X+A-E)*np.cos(alpha_e)
    y_e0 = M*(Z/2+X+A-E)*np.sin(alpha_e)
    # Parameter Range of Edge Round
    theta3_min = np.arctan((Y11[len(Y11)-1]-y_e0)/(X11[len(X11)-1]-x_e0))
    theta3_max = np.arctan((y_e-y_e0)/(x_e-x_e0))
    THETA3 = np.linspace(theta3_min,theta3_max,seg_edge_r) 
    X21 = M*E*np.cos(THETA3) + x_e0
    Y21 = M*E*np.sin(THETA3) + y_e0
    X22 = X21
    Y22 = -Y21
    X22 = X22[::-1]
    Y22 = Y22[::-1]

    # Root Round Curve of Tooth
    # Condition Check
    # alpha_ts = Start Angle of Root Round Curve
    # THETA_s = Substitution Variable to plot Root Round Curve
    alpha_ts = (2*(C*(1-np.sin(alpha_0))-D)*np.sin(alpha_0)+B)/(Z*np.cos(alpha_0)) - 2*C*np.cos(alpha_0)/Z + np.pi/(2*Z)
    theta_te = 2*C*np.cos(alpha_0)/Z - 2*(D-X-C*(1-np.sin(alpha_0)))*np.cos(alpha_0)/(Z*np.sin(alpha_0))
    THETA_t  = np.linspace(0,theta_te,seg_root_r)
    if (C!=0) and ((D-X-C)==0) :
        # mc를 반지름으로 하는 원호를 그려서 대체하게 됨
        THETA_s = (np.pi/2)*np.ones(len(THETA_t))
    elif (D-X-C)!=0 :
        THETA_s = np.arctan((M*Z*THETA_t/2)/(M*D-M*X-M*C))
    X31 = M*( (Z/2+X-D+C)*np.cos(THETA_t+alpha_ts) + (Z/2)*THETA_t*np.sin(THETA_t+alpha_ts) - C*np.cos(THETA_s+THETA_t+alpha_ts) )
    Y31 = M*( (Z/2+X-D+C)*np.sin(THETA_t+alpha_ts) - (Z/2)*THETA_t*np.cos(THETA_t+alpha_ts) - C*np.sin(THETA_s+THETA_t+alpha_ts) )
    X32 = X31
    Y32 = -Y31
    X32 = X32[::-1]
    Y32 = Y32[::-1]

    # Outer Arc
    THETA6 = np.linspace(alpha_e,alpha_m,SEG_OUTER) 
    X41 = M*(Z/2+A+X)*np.cos(THETA6)
    Y41 = M*(Z/2+A+X)*np.sin(THETA6)
    X42 = X41
    Y42 = -Y41
    X42 = X42[::-1]
    Y42 = Y42[::-1]

    # Root Arc
    THETA7 = np.linspace(0,alpha_ts,SEG_ROOT) 
    X51 = M*(Z/2-D+X)*np.cos(THETA7)
    Y51 = M*(Z/2-D+X)*np.sin(THETA7)
    X52 = X51
    Y52 = -Y51
    X52 = X52[::-1]
    Y52 = Y52[::-1]

    # Plot Whole Gear
    p_angle = 2*np.pi/z
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
    # Align to top
    align_angle = np.pi/2-np.pi/Z
    X2 = np.cos(align_angle)*X1 - np.sin(align_angle)*Y1
    Y2 = np.sin(align_angle)*X1 + np.cos(align_angle)*Y1
    # Array
    for i in range(0,Z) :
        # Rotate
        Xtemp = np.cos(p_angle*i)*X2 - np.sin(p_angle*i)*Y2
        Ytemp = np.sin(p_angle*i)*X2 + np.cos(p_angle*i)*Y2
        # Othogonal Transformation
        Xtemp = Xtemp + X_0
        Ytemp = Ytemp + Y_0
        # Plot
        plt.plot(Xtemp,Ytemp,'-',linewidth=1.5,color='black',label='_nolegend_')
    
    # Base Cicle
    THETA0 = np.linspace(0.0,2*np.pi,SEG_CIRCLE)
    base_dia = M*Z*np.cos(alpha_0)
    X_base = base_dia/2*np.sin(THETA0) + X_0
    Y_base = base_dia/2*np.cos(THETA0) + Y_0
    plt.plot(X_base, Y_base, ':', linewidth=1.0, color='cyan',label='Base Cicle')
    
    # Pitch Circle
    pitch_dia = M*Z
    X_pitch = pitch_dia/2*np.sin(THETA0) + X_0
    Y_pitch = pitch_dia/2*np.cos(THETA0) + Y_0
    plt.plot(X_pitch, Y_pitch, ':', linewidth=1.0, color='magenta',label='Pitch Circle')

    # Offset Circle
    offset_dia = 2*M*(Z/2+X)
    X_offset = (offset_dia/2)*np.sin(THETA0) + X_0
    Y_offset = (offset_dia/2)*np.cos(THETA0) + Y_0
    plt.plot(X_offset, Y_offset, ':', linewidth=1.0, color='red',label='Offset Circle')

    # Outer Circle
    outer_dia = 2*M*(Z/2+X+A)
    X_out = (outer_dia/2)*np.sin(THETA0) + X_0
    Y_out = (outer_dia/2)*np.cos(THETA0) + Y_0
    plt.plot(X_out, Y_out, ':', linewidth=1.0, color='brown',label='Outer Circle')

    # Root Circle
    root_dia = 2*M*(Z/2+X-D)
    X_root = (root_dia/2)*np.sin(THETA0) + X_0
    Y_root = (root_dia/2)*np.cos(THETA0) + Y_0
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
            
    # Save Figure
    #plt.legend(loc="lower left", ncol=1, fontsize="x-small")
    # Total Gear
    Result = os.path.join(WorkingDirectory, f'Result.png')
    plt.savefig(Result,dpi=100)
    #Result = os.path.join(WorkingDirectory, f'Result.svg')
    plt.savefig(Result)
    # One Tooth
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
    plt.savefig(Result2)

    ####################
    # DXF
    doc = ezdxf.new('R2000')
    msp = doc.modelspace()
    # Spline
    Xdxf1 = np.concatenate((X22,X13,X33))
    Ydxf1 = np.concatenate((Y22,Y13,Y33))
    Xdxf3 = np.cos(np.pi/2+np.pi/Z)*Xdxf1 - np.sin(np.pi/2+np.pi/Z)*Ydxf1
    Ydxf3 = np.sin(np.pi/2+np.pi/Z)*Xdxf1 + np.cos(np.pi/2+np.pi/Z)*Ydxf1
    Xdxf4 = -Xdxf3
    Ydxf4 =  Ydxf3
    Xdxf5 = Xdxf3 + X_0
    Ydxf5 = Ydxf3 + Y_0
    Xdxf6 = Xdxf4 + X_0
    Ydxf6 = Ydxf4 + Y_0
    cpoint5 = [(Xdxf5[0],Ydxf5[0])]
    cpoint6 = [(Xdxf6[0],Ydxf6[0])]
    for i in range(0,len(Xdxf5)) :
        cpoint5.append((Xdxf5[i],Ydxf5[i]))
        cpoint6.append((Xdxf6[i],Ydxf6[i]))
    msp.add_open_spline(cpoint5)
    msp.add_open_spline(cpoint6)
    
    # Arc outer
    outer_angle = np.arctan(Xdxf3[0]/Ydxf3[0])*360.0/(2.0*np.pi)
    msp.add_arc(center=(X_0,Y_0), radius=outer_dia/2, start_angle=90+outer_angle, end_angle=90-outer_angle)
    
    # Arc root
    root_start_angle = 180/Z
    root_end_angle = np.arctan(Xdxf3[len(Xdxf3)-1]/Ydxf3[len(Xdxf3)-1])*360.0/(2.0*np.pi)
    msp.add_arc(center=(X_0,Y_0), radius=root_dia/2, start_angle=90-root_start_angle, end_angle=90+root_end_angle)
    msp.add_arc(center=(X_0,Y_0), radius=root_dia/2, start_angle=90-root_end_angle, end_angle=90+root_start_angle)
    
    # Line dividing
    Xdivide1 = -np.sin(np.pi/Z)*root_dia/2 + X_0
    Ydivide1 =  np.cos(np.pi/Z)*root_dia/2 + Y_0
    Xdivide2 = -np.sin(-np.pi/Z)*root_dia/2 + X_0
    Ydivide2 =  np.cos(-np.pi/Z)*root_dia/2 + Y_0
    msp.add_line([X_0,Y_0],[Xdivide1,Ydivide1])
    msp.add_line([X_0,Y_0],[Xdivide2,Ydivide2])
    
    # Circle offset
    msp.add_circle((X_0,Y_0),radius=offset_dia/2)

    Result = os.path.join(values['-WorkingDirectoty-'], f'Result.dxf')
    doc.saveas(Result)
    
    ####################
    # Write Gear Spec
    Result3 = os.path.join(values['-WorkingDirectoty-'], f'Result.csv')
    fileout = open(Result3, "w")
    if alpha == 20 :
        fileout.write("Type,"+"Standard"+",\n")
    else :
        fileout.write("Type,"+"Non-Standard"+",\n")
    fileout.write("Module,"+repr(M)+",mm"+"\n")
    fileout.write("Pressure Angle,"+repr(ALPHA)+",deg"+"\n")
    fileout.write("Teeth Number,"+repr(Z)+",ea"+"\n")
    fileout.write("Offset Factor,"+repr(X)+",\n")
    fileout.write("Offset,"+repr(X*M)+",mm"+"\n")
    fileout.write("Backlash Factor,"+repr(B)+",\n")
    fileout.write("Backlash,"+repr(B*M)+",mm"+"\n")
    fileout.write("Addendum Factor,"+repr(A)+","+"\n")
    fileout.write("Addendum,"+repr(A*M)+",mm"+"\n")
    fileout.write("Dedendum Factor,"+repr(D)+","+"\n")
    fileout.write("Dedendum,"+repr(D*M)+",mm"+"\n")
    fileout.write("Total Tooth Height,"+repr(A+D)+",mm"+"\n")
    fileout.write("Base Circle Dia,"+repr(base_dia)+",mm"+"\n")
    fileout.write("Pitch Circle Dia,"+repr(pitch_dia)+",mm"+"\n")
    fileout.write("Offset Circle Dia,"+repr(offset_dia)+",mm"+"\n")
    fileout.write("Root Circle Dia,"+repr(root_dia)+",mm"+"\n")
    fileout.write("Outer Circle Dia,"+repr(outer_dia)+",mm"+"\n")
    fileout.close()    

##############################
# GUI
sg.theme('Default')

left_col = [[sg.Text('1. Gear Spec',font='ARIAL 16')],
            [sg.Text('Module, m =',size = (32,1)),sg.Input(1.0,key='-m-',size = (10,1)),sg.Text('[mm], (>0)')],
            [sg.Text('Teeth Number, z =',size = (32,1)),sg.Input(15,key='-z-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Pressure Angle [Deg], alpha =',size = (32,1)),sg.Input(20.0,key='-alpha-',size = (10,1)),sg.Text('[deg]')],
            [sg.Text('Offset Factor, x =',size = (32,1)),sg.Input(0.2,key='-x-',size = (10,1)),sg.Text('(-1~+1)')],
            [sg.Text('Backlash Factor, b =',size = (32,1)),sg.Input(0.05,key='-b-',size = (10,1)),sg.Text('(0~1)')],
            [sg.Text('Addendum Factor, a =',size = (32,1)),sg.Input(1.0,key='-a-',size = (10,1)),sg.Text('(0~1)')],
            [sg.Text('Dedendum Factor, d =',size = (32,1)),sg.Input(1.25,key='-d-',size = (10,1)),sg.Text('(0~1)')],
            [sg.Text('Radius Factor of Edge Round of Hob, c =',size = (32,1)),sg.Input(0.2,key='-c-',size = (10,1))],
            [sg.Text('Radius Factor of Edge Round of Tooth, e =',size = (32,1)),sg.Input(0.2,key='-e-',size = (10,1))],

            [sg.Text('2. Graphics',font='ARIAL 16')],
            [sg.Text('Center of Gear, x_0 =',size = (32,1)),sg.Input(0.0,key='-x_0-',size = (10,1)),sg.Text('[mm]')],
            [sg.Text('Center of Gear, y_0 =',size = (32,1)),sg.Input(0.0,key='-y_0-',size = (10,1)),sg.Text('[mm]')],
            [sg.Text('Segmentation Numbers, seg_circle =',size = (32,1)),sg.Input(360,key='-seg_circle-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Segmentation Numbers, seg_involute =',size = (32,1)),sg.Input(15,key='-seg_involute-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Segmentation Numbers, seg_edge_r =',size = (32,1)),sg.Input(5,key='-seg_edge_r-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Segmentation Numbers, seg_root_r =',size = (32,1)),sg.Input(5,key='-seg_root_r-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Segmentation Numbers, seg_center =',size = (32,1)),sg.Input(5,key='-seg_center-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Segmentation Numbers, seg_outer =',size = (32,1)),sg.Input(5,key='-seg_outer-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Segmentation Numbers, seg_root =',size = (32,1)),sg.Input(5,key='-seg_root-',size = (10,1)),sg.Text('[ea]')],
            [sg.Text('Scale for One Tooth, scale =',size = (32,1)),sg.Input(0.7,key='-scale-',size = (10,1)),sg.Text('(0.1~1)')]]

right_col = [[sg.Text('Working Directory :',size=(8,1)), sg.Input('./Result/',key='-WorkingDirectoty-',size=(16,1)), sg.FolderBrowse()],
            [sg.Image('FGPG2.png',size=(500,500),key='-IMAGE-')],
            [sg.Text('Hello',key='-TEXT-')],
            [sg.Button('Load'), sg.Button('Run'), sg.Button('Toggle'), sg.Button('Exit')]]

CurrentImage = 'Result'

layout = [[sg.Column(left_col), sg.VSeperator(), sg.Column(right_col)]]

window = sg.Window('FGPG2', layout)

while True:
    event, values = window.read()

    try:
        m = float(values['-m-'])
        z = int(values['-z-'])
        alpha = float(values['-alpha-'])
        x = float(values['-x-'])
        b = float(values['-b-'])
        a = float(values['-a-'])
        d = float(values['-d-'])
        c = float(values['-c-'])
        e = float(values['-e-'])
        x_0 = float(values['-x_0-'])
        y_0 = float(values['-y_0-'])
        seg_circle = int(values['-seg_circle-'])
        seg_involute = int(values['-seg_involute-'])
        seg_edge_r = int(values['-seg_edge_r-'])
        seg_root_r = int(values['-seg_root_r-'])
        seg_center = int(values['-seg_center-'])
        seg_outer = int(values['-seg_outer-'])
        seg_root = int(values['-seg_root-'])
        scale = float(values['-scale-'])
        WorkingDirectory = values['-WorkingDirectoty-']
    except:
        sg.popup('Type error.')

    if event in (sg.WIN_CLOSED, 'Exit'):
        break
    elif event == 'Load':
        Inputs = os.path.join(WorkingDirectory, f'Inputs.dat')
        if os.path.exists(Inputs) :
            window.load_from_disk(Inputs)
            window['-TEXT-'].update('Load : OK')
        else :
            sg.popup('The File not exists : %s'%Inputs)
    elif event == 'Run':
        if os.path.exists(WorkingDirectory) :
            FGPG2_PLOT(m,z,alpha,x,b,a,d,c,e,x_0,y_0,seg_circle,seg_involute,seg_edge_r,seg_root_r,seg_center,seg_outer,seg_root,scale)
            Result = os.path.join(WorkingDirectory, f'Result.png')
            window['-IMAGE-'].update(Result,size=(500,500))
            Inputs = os.path.join(WorkingDirectory, f'Inputs.dat')
            window.save_to_disk(Inputs)
            window['-TEXT-'].update('Run : OK')
        else :
            sg.popup('The Directory not exists : %s'%WorkingDirectory)
    elif event == 'Toggle':
        if CurrentImage == 'Result':
            Result = os.path.join(WorkingDirectory, f'Result2.png')
            CurrentImage = 'Result2'
            window['-IMAGE-'].update(Result,size=(500,500))
            window['-TEXT-'].update('Toggle to One Tooth : OK')
        elif CurrentImage == 'Result2':
            Result = os.path.join(WorkingDirectory, f'Result.png')
            CurrentImage = 'Result'
            window['-IMAGE-'].update(Result,size=(500,500))
            window['-TEXT-'].update('Toggle to Whole Gear : OK')

window.close()
