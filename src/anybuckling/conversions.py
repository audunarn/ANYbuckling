'''
Stress/force conversions for cylindrical and conical shells.

Transferred from ANYstructure (anystruct/helper.py).
'''
import math

def helper_cylinder_stress_to_force_to_stress(stresses = None, forces = None, geometry = None, shell_t = 0,
                                              shell_radius = 0, shell_spacing = 0,
                                              hw = 0, tw = 0, b = 0, tf = 0, CylinderAndCurvedPlate = None,
                                              conical = False, psd = 0, cone_r1 = 0, cone_r2 = 0, cone_alpha = 0,
                                              shell_lenght_l = 0):

    hw = 0 if hw is None else hw
    tw = 0 if tw is None else tw
    b = 0 if b is None else b
    tf = 0 if tf is None else tf
    has_longitudinal_stiffener = geometry in [3, 4, 7, 8] and shell_spacing not in [None, 0] and \
        hw * tw + b * tf > 0

    A = hw * tw + b * tf if has_longitudinal_stiffener else 0
    eq_thk = shell_t + A/shell_spacing if has_longitudinal_stiffener else shell_t

    Itot = CylinderAndCurvedPlate.get_Itot(hw=hw if has_longitudinal_stiffener else 0,
                                           tw=tw if has_longitudinal_stiffener else 0,
                                           b=b if has_longitudinal_stiffener else 0,
                                           tf=tf if has_longitudinal_stiffener else 0,
                                           r=shell_radius,
                                           s=shell_spacing,
                                           t=shell_t)

    if forces is not None and stresses is None:
        if not conical:
            Nsd, Msd, Tsd, Qsd = forces
            sasd = (Nsd / 2) / (math.pi * shell_radius * eq_thk) * 1000
            smsd = (Msd/ Itot) * \
                   (shell_radius + shell_t / 2) * 1000000
            tTsd = (Tsd* 10 ** 6) / (2 * math.pi * shell_t * math.pow(shell_radius, 2))
            tQsd = Qsd / (math.pi * shell_radius * shell_t) * 1000
            shsd = 0
            return sasd, smsd, tTsd, tQsd, shsd
        else:
            Nsd, M1sd, M2sd, Tsd, Q1sd, Q2sd = forces
            r = shell_radius if shell_radius not in [None, 0] else (cone_r1 + cone_r2) / 2
            cos_alpha = math.cos(math.radians(cone_alpha))
            te = shell_t * cos_alpha
            pressure_factor = 1e6 if max(abs(shell_t), abs(r), abs(shell_lenght_l)) < 100 else 1
            pressure = psd * pressure_factor
            sasd = pressure * r / (2 * te) + Nsd / (2 * math.pi * r * te) * 1000
            smsd = math.sqrt(math.pow(M1sd, 2) + math.pow(M2sd, 2)) * 1000 / (math.pi * math.pow(r, 2) * te)
            shsd = pressure * r / te
            tTsd = Tsd * 1000 / (2 * math.pi * math.pow(r, 2) * te)
            tQsd = math.sqrt(math.pow(Q1sd, 2) + math.pow(Q2sd, 2)) * 1000 / (math.pi * r * te)
            return sasd, smsd, tTsd, tQsd, shsd

    else:
        if not conical:
            sasd, smsd, tTsd, tQsd, shsd = stresses
            Nsd = (sasd * 2 * math.pi * shell_radius * eq_thk) / 1000
            Msd = (smsd * Itot) / ((shell_radius + shell_t / 2) * 1000000)
            Tsd = tTsd * 2 * math.pi * shell_t * math.pow(shell_radius, 2) / 1000000
            Qsd = tQsd * math.pi * shell_radius * shell_t / 1000
        else:
            r = shell_radius if shell_radius not in [None, 0] else (cone_r1 + cone_r2) / 2
            cos_alpha = math.cos(math.radians(cone_alpha))
            te = shell_t * cos_alpha
            pressure_factor = 1e6 if max(abs(shell_t), abs(r), abs(shell_lenght_l)) < 100 else 1
            pressure = psd * pressure_factor
            sasd, smsd, tTsd, tQsd, shsd = stresses
            Nsd = (sasd - pressure * r / (2 * te)) * (2 * math.pi * r * te) / 1000
            M1sd = smsd * math.pi * math.pow(r, 2) * te / 1000
            M2sd = 0
            Tsd = tTsd * 2 * math.pi * math.pow(r, 2) * te / 1000
            Q1sd = tQsd * math.pi * r * te / 1000
            Q2sd = 0
            return Nsd, M1sd, M2sd, Tsd, Q1sd, Q2sd, shsd

        return Nsd, Msd, Tsd, Qsd, shsd


if __name__ == '__main__':
    from tkinter import *

    class AllTkinterWidgets:
        def __init__(self, master):
            frame = Frame(master, width=500, height=400, bd=1)
            frame.pack()

            iframe5 = Frame(frame, bd=2, relief=RAISED)
            iframe5.pack(expand=1, fill=X, pady=10, padx=5)
            c = Canvas(iframe5, bg='white', width=340, height=200)
            c.pack()

            height = 150
            radius = 150
            offset_oval = 30
            start_x_cyl = 150
            start_y_cyl = 20
            coord1 = start_x_cyl, start_y_cyl, start_x_cyl + radius, offset_oval
            coord2 = start_x_cyl, start_y_cyl + height, start_x_cyl + radius, offset_oval+ height

            arc_1 = c.create_oval(coord1, width = 5, fill = 'grey90')
            arc_2 = c.create_arc(coord2, extent = 180, start = 180,style=ARC, width = 3)

            line1 = c.create_line(coord1[0], coord1[1]+offset_oval/4,
                                  coord1[0], coord1[1]+height+offset_oval/4,
                                  width = 3)
            line2 = c.create_line(coord1[0]+radius, coord1[1]+offset_oval/4,
                                  coord1[0]+radius, coord1[1]+height+offset_oval/4,
                                  width = 3)
            num_stf = 10
            for line_num in range(1,num_stf,1):
                angle = 180 - 180/(num_stf) *line_num
                arc_x, arc_y = 1*math.cos(math.radians(angle)), 0.5*math.sin(math.radians(angle))
                arc_x = (arc_x + 1)/2

                line1 = c.create_line(coord1[0] + radius*arc_x,
                                      coord1[1] +2*arc_y*offset_oval/3,
                                      coord1[0] + radius*arc_x,
                                      coord1[1] + height +2*arc_y*offset_oval/3,fill = 'blue')
            num_ring_stiff = 5
            for ring_stf in range(1,num_ring_stiff+1,1):
                coord3 = coord1[0], coord1[1]+(height/(num_ring_stiff+1))*ring_stf,  \
                         start_x_cyl +radius, coord1[3]+ (height/(num_ring_stiff+1))*ring_stf,
                arc_2 = c.create_arc(coord3, extent=180, start=180, style=ARC, width=2,fill = 'orange', outline = 'orange')

            num_ring_girder = 1
            for ring_girder in range(1, num_ring_girder+1,1):
                coord3 = coord1[0], coord1[1]+(height/(num_ring_girder+1))*ring_girder,  \
                         start_x_cyl+ radius, coord1[3]+ (height/(num_ring_girder+1))*ring_girder,
                arc_2 = c.create_arc(coord3, extent=180, start=180, style=ARC, width=4, fill = 'grey', outline = 'grey')

            iframe5.pack(expand=1, fill=X, pady=10, padx=5)


    root = Tk()
    # root.option_add('*font', ('verdana', 10, 'bold'))
    all = AllTkinterWidgets(root)
    root.title('Tkinter Widgets')
    root.mainloop()
