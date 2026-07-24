'''
Prescriptive buckling and scantling checks for flat stiffened panels
according to DNV standards (DNV-RP-C201 / DNV-OS-C101).

Transferred from ANYstructure (anystruct/calc_structure.py).
Contains the panel data model (Structure), scantling checks
(CalcScantlings) and the stiffened flat-plate buckling calculation
(AllStructure).
'''
import math
import numpy as np
from scipy.optimize import minimize

class Structure():
    '''
    Setting the properties for the plate and the stiffener. Takes a dictionary as argument.
    '''
    def __init__(self, main_dict: dict = None, *args, **kwargs):
        super(Structure,self).__init__()
        if main_dict is None:
            self._panel_or_shell = None
            self._plate_th = None
            self._web_height = None
            self._web_th = None
            self._flange_width = None
            self._flange_th = None
            self._mat_yield = None
            self._mat_factor = None
            self._span = None
            self._spacing = None
            self._structure_type = None
            self._sigma_y1 = None
            self._sigma_y2 = None
            self._sigma_x1 = None
            self._sigma_x2 = None
            self._tauxy = None
            self._plate_kpp = None
            self._stf_kps = None
            self._km1 = None
            self._km2 = None
            self._km3 = None
            self._stiffener_type = None
            self._structure_types = None
            self._dynamic_variable_orientation = None
            if self._structure_type is None:
                self._dynamic_variable_orientation = None
            self._puls_method = None
            self._puls_boundary = None
            self._puls_stf_end = None
            self._puls_sp_or_up = None
            self._puls_up_boundary = None
            self._zstar_optimization = None
            try:
                self._girder_lg = None
            except KeyError:
                self._girder_lg = None
            try:
                self._pressure_side = None
            except KeyError:
                self._pressure_side = None
            self._panel_or_shell = None
        else:
            self._main_dict = main_dict
            if 'panel or shell' not in main_dict.keys():
                self._panel_or_shell = 'panel'
            else:
                self._panel_or_shell = main_dict['panel or shell'][0]
            self._plate_th = main_dict['plate_thk'][0]
            self._web_height = main_dict['stf_web_height'][0]
            self._web_th = main_dict['stf_web_thk'][0]
            self._flange_width = main_dict['stf_flange_width'][0]
            self._flange_th = main_dict['stf_flange_thk'][0]
            self._mat_yield = main_dict['mat_yield'][0]
            self._mat_factor = main_dict['mat_factor'][0]
            self._span = main_dict['span'][0]
            self._spacing = main_dict['spacing'][0]
            self._structure_type = main_dict['structure_type'][0]
            self._sigma_y1=main_dict['sigma_y1'][0]
            self._sigma_y2=main_dict['sigma_y2'][0]
            self._sigma_x1 = main_dict['sigma_x1'][0]
            self._sigma_x2 = main_dict['sigma_x2'][0]
            self._tauxy=main_dict['tau_xy'][0]
            self._plate_kpp = main_dict['plate_kpp'][0]
            self._stf_kps = main_dict['stf_kps'][0]
            self._km1 = main_dict['stf_km1'][0]
            self._km2 = main_dict['stf_km2'][0]
            self._km3 = main_dict['stf_km3'][0]
            self._stiffener_type=main_dict['stf_type'][0]
            self._structure_types = main_dict['structure_types'][0]
            self._dynamic_variable_orientation = None
            if self._structure_type in self._structure_types['vertical']:
                self._dynamic_variable_orientation = 'z - vertical'
            elif self._structure_type in self._structure_types['horizontal']:
                self._dynamic_variable_orientation = 'x - horizontal'
            self._puls_method = main_dict['puls buckling method'][0]
            self._puls_boundary = main_dict['puls boundary'][0]
            self._puls_stf_end = main_dict['puls stiffener end'][0]
            self._puls_sp_or_up = main_dict['puls sp or up'][0]
            self._puls_up_boundary = main_dict['puls up boundary'][0]

            self._zstar_optimization = main_dict['zstar_optimization'][0]
            try:
                self._girder_lg=main_dict['girder_lg'][0]
            except KeyError:
                self._girder_lg = 10
            try:
                self._pressure_side = main_dict['press_side'][0]
            except KeyError:
                self._pressure_side = 'both sides'
            self._panel_or_shell = main_dict['panel or shell'][0]

    # Property decorators are used in buckling of shells. IN mm!
    @property # in mm
    def hw(self):
        assert self._web_height is not None, 'Variable missing: self._web_height'
        return self._web_height * 1000
    @hw.setter # in mm
    def hw(self, val):
        self._web_height = val / 1000
    @property # in mm
    def tw(self):
        assert self._web_th is not None, 'Variable missing: self._web_th - web thickness'
        return self._web_th * 1000
    @tw.setter # in mm
    def tw(self, val):
        self._web_th = val / 1000
    @property # in mm
    def b(self):
        assert self._flange_width is not None, 'Variable missing: self._flange_width'
        return self._flange_width * 1000
    @b.setter # in mm
    def b(self, val):
        self._flange_width = val / 1000
    @property # in mm
    def tf(self):
        assert self._flange_th is not None, 'Variable missing: self._flange_th'
        return self._flange_th * 1000
    @tf.setter # in mm
    def tf(self, val):
        self._flange_th = val / 1000
    @property  # in mm
    def spacing(self):
        assert self._spacing is not None, 'Variable missing: self._spacing'
        return self._spacing* 1000
    @spacing.setter  # in mm
    def spacing(self, val):
        self._spacing = val / 1000
    @property  # in mm
    def t(self):
        assert self._plate_th is not None, 'Variable missing: self._plate_th'
        return self._plate_th* 1000
    @t.setter  # in mm
    def t(self, val):
        self._plate_th = val / 1000
    @property  # in mm
    def panel_or_shell(self):

        return self._panel_or_shell
    @panel_or_shell.setter  # in mm
    def panel_or_shell(self, val):
        self._panel_or_shell = val   
    @property
    def stiffener_type(self):
        return self._stiffener_type
    @stiffener_type.setter
    def stiffener_type(self, val):
        self._stiffener_type = val

    @property
    def E(self):
        return self._E
    @E.setter
    def E(self, val):
        self._E = val
    @property
    def v(self):
        return self._v
    @v.setter
    def v(self, val):
        self._v = val
    @property
    def mat_yield(self):
        return self._mat_yield
    @mat_yield.setter
    def mat_yield(self, val):
        self._mat_yield = val

    @property
    def mat_factor(self):
        assert self._mat_factor is not None, 'Missing variable: self._mat_factor'
        return self._mat_factor
    @mat_factor.setter
    def mat_factor(self, val):
        self._mat_factor = val

    @property
    def span(self):
        assert self._span is not None, 'Missing variable: self._span: span of stiffener'
        return self._span
    @span.setter
    def span(self, val):
        self._span = val
        
    @property
    def sigma_x1(self):
        assert self._sigma_x1 is not None, 'Missing variable: self._sigma_x1'
        return self._sigma_x1
    @sigma_x1.setter
    def sigma_x1(self, val):
        self._sigma_x1 = val
    @property
    def sigma_x2(self):
        assert self._sigma_x2 is not None, 'Missing variable: self._sigma_x2'
        return self._sigma_x2
    @sigma_x2.setter
    def sigma_x2(self, val):
        self._sigma_x2 = val
        
    @property
    def sigma_y1(self):
        assert self._sigma_y1 is not None, 'Missing variable: self._sigma_y1'
        return self._sigma_y1
    @sigma_y1.setter
    def sigma_y1(self, val):
        self._sigma_y1 = val
    
    @property
    def sigma_y2(self):
        assert self._sigma_y2 is not None, 'Missing variable: self._sigma_y2'
        return self._sigma_y2
    @sigma_y2.setter
    def sigma_y2(self, val):
        self._sigma_y2 = val
        
    @property
    def tau_xy(self):
        assert self._tauxy is not None, 'Missing variable: self._tauxy'
        return self._tauxy
    @tau_xy.setter
    def tau_xy(self, val):
        self._tauxy = val
        
    @property
    def girder_lg(self):
        assert self._girder_lg is not None, 'Missing variable: self._girder_lg'
        return self._girder_lg
    @girder_lg.setter
    def girder_lg(self, val):
        self._girder_lg = val


    def __str__(self):
        '''
        Returning all properties.
        '''
        try:
            return \
                str(
                '\n Plate field span:              ' + str(round(self._span*1000)) + ' mm' +
                '\n Stiffener spacing:             ' + str(self._spacing*1000)+' mm'+
                '\n Plate thickness:               ' + str(self._plate_th*1000)+' mm'+
                '\n Stiffener web height:          ' + str(self._web_height*1000)+' mm'+
                '\n Stiffener web thickness:       ' + str(self._web_th*1000)+' mm'+
                '\n Stiffener flange width:        ' + str(self._flange_width*1000)+' mm'+
                '\n Stiffener flange thickness:    ' + str(self._flange_th*1000)+' mm'+
                '\n Material yield:                ' + str(self._mat_yield/1e6)+' MPa'+
                '\n Structure/stiffener type:      ' + str(self._structure_type)+'/'+(self._stiffener_type)+
                '\n Dynamic load varible_          ' + str(self._dynamic_variable_orientation)+
                '\n Plate fixation paramter,kpp:   ' + str(self._plate_kpp) + ' ' +
                '\n Stf. fixation paramter,kps:    ' + str(self._stf_kps) + ' ' +
                '\n Global stress, sig_y1/sig_y2:  ' + str(round(self._sigma_y1,3))+'/'+str(round(self._sigma_y2,3))+ ' MPa' +
                '\n Global stress, sig_x1/sig_x2:   ' + str(round(self._sigma_x1,3))+'/'+str(round(self._sigma_x2,3))+ ' MPa' +
                '\n Global shear, tau_xy:          ' + str(round(self._tauxy,3)) + ' MPa' +
                '\n km1,km2,km3:                   ' + str(self._km1)+'/'+str(self._km2)+'/'+str(self._km3)+
                '\n Pressure side (p-plate/s-stf): ' + str(self._pressure_side) + ' ')
        except TypeError:
            return \
                str(
                    '\n Stiffener spacing:             ' + str(self._spacing * 1000) + ' mm' +
                    '\n Plate thickness:               ' + str(self._plate_th * 1000) + ' mm' +
                    '\n Stiffener web height:          ' + str(self._web_height * 1000) + ' mm' +
                    '\n Stiffener web thickness:       ' + str(self._web_th * 1000) + ' mm' +
                    '\n Stiffener flange width:        ' + str(self._flange_width * 1000) + ' mm' +
                    '\n Stiffener flange thickness:    ' + str(self._flange_th * 1000) + ' mm' )

    def get_beam_string(self, short = False):
        ''' Returning a string. '''
        if type(self._stiffener_type) != str:
            print('error')

        base_name = self._stiffener_type+ '_' + str(round(self._web_height*1000, 0)) + 'x' + \
                   str(round(self._web_th*1000, 0))
        if self._stiffener_type == 'FB':
            ret_str = base_name
        elif self._stiffener_type in ['L-bulb', 'bulb', 'hp']:
            if not short:
                ret_str = 'Bulb'+str(int(self._web_height*1000 + self._flange_th*1000))+'x'+\
                          str(round(self._web_th*1000, 0))+ '_(' +str(round(self._web_height*1000, 0)) + 'x' + \
                       str(round(self._web_th*1000, 0))+'_'+ str(round(self._flange_width*1000, 0)) + 'x' + \
                          str(round(self._flange_th*1000, 0))+')'
            else:
                ret_str = 'Bulb'+str(int(self._web_height*1000 + self._flange_th*1000))+'x'+\
                      str(round(self._web_th*1000, 0))
        else:
            ret_str = base_name + '__' + str(round(self._flange_width*1000, 0)) + 'x' + \
                      str(round(self._flange_th*1000, 0))

        ret_str = ret_str.replace('.', '_')

        return ret_str
        # base_name = self._stiffener_type+ '_' + str(round(self._web_height*1000, 0)) + 'x' + \
        #            str(round(self._web_th*1000, 0))
        # if self._stiffener_type == 'FB':
        #     ret_str = base_name
        # else:
        #     ret_str = base_name + '__' + str(round(self._flange_width*1000, 0)) + 'x' + \
        #               str(round(self._flange_th*1000, 0))
        #
        # ret_str = ret_str.replace('.', '_')
        #
        # return ret_str

    def get_structure_types(self):
        return self._structure_types

    def get_z_opt(self):
        return self._zstar_optimization

    def get_puls_method(self):
        return self._puls_method

    def get_puls_boundary(self):
        return self._puls_boundary

    def get_puls_stf_end(self):
        return self._puls_stf_end

    def get_puls_sp_or_up(self):
        return self._puls_sp_or_up

    def get_puls_up_boundary(self):
        return self._puls_up_boundary

    def get_one_line_string(self):
        ''' Returning a one line string. '''
        return 'pl_'+str(round(self._spacing*1000, 1))+'x'+str(round(self._plate_th*1000,1))+' stf_'+self._stiffener_type+\
               str(round(self._web_height*1000,1))+'x'+str(round(self._web_th*1000,1))+'+'\
               +str(round(self._flange_width*1000,1))+'x'+\
               str(round(self._flange_th*1000,1))

    def get_report_stresses(self):
        'Return the stresses to the report'
        return 'sigma_y1: '+str(round(self._sigma_y1,1))+' sigma_y2: '+str(round(self._sigma_y2,1))+ \
               ' sigma_x1: ' + str(round(self._sigma_x1,1)) +' sigma_x2: ' + str(round(self._sigma_x2,1))+\
               ' tauxy: '+ str(round(self._tauxy,1))

    def get_extended_string(self):
        ''' Some more information returned. '''
        return 'span: '+str(round(self._span,4))+' structure type: '+ self._structure_type + ' stf. type: ' + \
               self._stiffener_type + ' pressure side: ' + self._pressure_side

    def get_s(self):
        '''
        Return the spacing
        :return:
        '''
        return self._spacing
    def get_pl_thk(self):
        '''
        Return the plate thickness
        :return:
        '''
        return self._plate_th
    def get_web_h(self):
        '''
        Return the web heigh
        :return:
        '''
        return self._web_height
    def get_web_thk(self):
        '''
        Return the spacing
        :return:
        '''
        return self._web_th
    def get_fl_w(self):
        '''
        Return the flange width
        :return:
        '''
        return self._flange_width
    def get_fl_thk(self):
        '''
        Return the flange thickness
        :return:
        '''
        return self._flange_th
    def get_fy(self):
        '''
        Return material yield
        :return:
        '''
        return self._mat_yield

    def get_span(self):
        '''
        Return the span
        :return:
        '''
        return self._span
    def get_kpp(self):
        '''
        Return var
        :return:
        '''
        return self._plate_kpp
    def get_kps(self):
        '''
        Return var
        :return:
        '''
        return self._stf_kps
    def get_km1(self):
        '''
        Return var
        :return:
        '''
        return self._km1
    def get_km2(self):
        '''
        Return var
        :return:
        '''
        return self._km2
    def get_km3(self):
        '''
        Return var
        :return:
        '''
        return self._km3
    def get_side(self):
        '''
        Return the checked pressure side.
        :return: 
        '''
        return self._pressure_side
    def get_tuple(self):
        ''' Return a tuple of the plate stiffener'''
        return (self._spacing, self._plate_th, self._web_height, self._web_th, self._flange_width,
                self._flange_th, self._span, self._girder_lg, self._stiffener_type)

    def get_section_modulus(self, efficient_se = None, dnv_table = False):
        '''
        Returns the section modulus.
        :param efficient_se: 
        :return: 
        '''
        #Plate. When using DNV table, default values are used for the plate
        b1 = self._spacing if efficient_se==None else efficient_se
        tf1 = self._plate_th

        #Stiffener
        tf2 = self._flange_th
        b2 = self._flange_width
        h = self._flange_th+self._web_height+self._plate_th
        tw = self._web_th
        hw = self._web_height

        # cross section area
        Ax = tf1 * b1 + tf2 * b2 + hw * tw

        assert Ax != 0, 'Ax cannot be 0'
        # distance to center of gravity in z-direction
        ez = (tf1 * b1 * tf1 / 2 + hw * tw * (tf1 + hw / 2) + tf2 * b2 * (tf1 + hw + tf2 / 2)) / Ax

        #ez = (tf1 * b1 * (h - tf1 / 2) + hw * tw * (tf2 + hw / 2) + tf2 * b2 * (tf2 / 2)) / Ax
        # moment of inertia in y-direction (c is centroid)

        Iyc = (1 / 12) * (b1 * math.pow(tf1, 3) + b2 * math.pow(tf2, 3) + tw * math.pow(hw, 3))
        Iy = Iyc + (tf1 * b1 * math.pow(tf1 / 2, 2) + tw * hw * math.pow(tf1+hw / 2, 2) +
             tf2 * b2 * math.pow(tf1+hw+tf2 / 2, 2)) - Ax * math.pow(ez, 2)

        # elastic section moduluses y-axis
        Wey1 = Iy / (h - ez)
        Wey2 = Iy / ez

        return Wey1, Wey2
    def get_plasic_section_modulus(self):
        '''
        Returns the plastic section modulus
        :return:
        '''
        tf1 = self._plate_th
        tf2 = self._flange_th
        b1 = self._spacing
        b2 = self._flange_width
        h = self._flange_th+self._web_height+self._plate_th
        tw = self._web_th
        hw = self._web_height

        Ax = tf1 * b1 + tf2 * b2 + (h-tf1-tf2) * tw

        ezpl = (Ax/2-b1*tf1)/tw+tf1

        az1 = h-ezpl-tf1
        az2 = ezpl-tf2

        Wy1 = b1*tf1*(az1+tf1/2) + (tw/2)*math.pow(az1,2)
        Wy2 = b2*tf2*(az2+tf2/2)+(tw/2)*math.pow(az2,2)

        return Wy1+Wy2
    def get_shear_center(self):
        '''
        Returning the shear center
        :return:
        '''
        tf1 = self._plate_th
        tf2 = self._flange_th
        b1 = self._spacing
        b2 = self._flange_width
        h = self._flange_th+self._web_height+self._plate_th
        tw = self._web_th
        hw = self._web_height
        Ax = tf1 * b1 + tf2 * b2 + (h-tf1-tf2) * tw
        # distance to center of gravity in z-direction
        ez = (b2*tf2*tf2/2 + tw*hw*(tf2+hw/2)+tf1*b1*(tf2+hw+tf1/2)) / Ax

        # Shear center:
        # moment of inertia, z-axis
        Iz1 = tf1 * math.pow(b1, 3)
        Iz2 = tf2 * math.pow(b2, 3)
        ht = h - tf1 / 2 - tf2 / 2
        return (Iz1 * ht) / (Iz1 + Iz2) + tf2 / 2 - ez

    def get_moment_of_intertia_hp(self):
        return self.get_moment_of_intertia()

    def get_moment_of_intertia(self, efficent_se=None, only_stf = False, tf1 = None, reduced_tw = None,
                               plate_thk = None, plate_spacing = None):
        '''
        Returning moment of intertia.
        :return:
        '''

        if only_stf:
            tf1 = t = 0
            b1 = s_e = 0
        else:
            # if tf1 == None:
            #     tf1 = self._plate_th
            tf1 = t =  self._plate_th if tf1 == None else tf1
            b1 = s_e =self._spacing if efficent_se==None else efficent_se

        e_f = 0

        h = self._flange_th+self._web_height+tf1
        tw = self._web_th if reduced_tw == None else reduced_tw/1000
        hw = self._web_height
        tf2 = tf = self._flange_th
        b2 = bf = self._flange_width

        Ax = tf1 * b1 + tf2 * b2 + hw * tw
        Iyc = (1 / 12) * (b1 * math.pow(tf1, 3) + b2 * math.pow(tf2, 3) + tw * math.pow(hw, 3))
        ez = (tf1 * b1 * (h - tf1 / 2) + hw * tw * (tf2 + hw / 2) + tf2 * b2 * (tf2 / 2)) / Ax
        Iy = Iyc + (tf1 * b1 * math.pow(tf2 + hw + tf1 / 2, 2) + tw * hw * math.pow(tf2 + hw / 2, 2) +
             tf2 * b2 * math.pow(tf2 / 2, 2)) - Ax * math.pow(ez, 2)

        # ###
        # z_c = bf * tf * e_f / (s_e * t + hw * tw + bf * tf)
        # I_z = 1.0 / 12.0 * t * math.pow(s_e,3) + 1.0 / 12.0 * hw * math.pow(tw,3) + 1.0 / 12.0 * tf * math.pow(bf,3) +\
        #       t * s_e * math.pow(z_c,2) + \
        #       tw * hw * math.pow(z_c,2) + bf * tf * math.pow(e_f - z_c,2)
        # ###
        #
        # z_c = (bf * tf * (tf / 2.0 + t / 2.0 + hw) + hw * tw * (hw / 2.0 + t / 2.0)) / (s_e * t + hw * tw + bf * tf)
        # I_sef = 1.0 / 12.0 * tw * hw ** 3 + 1.0 / 12.0 * bf * tf ** 3 + 1.0 / 12.0 * s_e * t ** 3 + tw * hw * (
        #             hw / 2.0 + t / 2.0 - z_c) ** 2 + tf * bf * (hw + t / 2.0 + tf / 2.0 - z_c) ** 2 + s_e * t * z_c ** 2
        # print(I_sef, I_z, Iy)
        #print(2*(bf*(h/2)**3/12 + tf*(h-tf)**3/12) + tw*h**3/12, Iy)
        return Iy

    def get_Iz_moment_of_inertia(self, reduced_tw = None):
        tw = self._web_th*1000 if reduced_tw is None else reduced_tw
        hw = self._web_height * 1000
        tf2 = self._flange_th * 1000
        b2 = self._flange_width * 1000

        if self._stiffener_type == 'FB':
            Iz = math.pow(tw,3)*hw/12
        elif self._stiffener_type == 'T':
            Iz = hw*math.pow(tw,3)/12 + tf2*math.pow(b2,3)/12
        else:
            Czver = tw/2
            Czhor = b2/2
            Aver = hw*tw
            Ahor = b2*tf2
            Atot = Aver+Ahor

            Czoverall = Aver*Czver/Atot + Ahor*Czhor/Atot
            dz = Czver - Czoverall

            Iver = (1/12)*hw*math.pow(tw,3) + Aver*math.pow(dz,2)

            dz = Czhor-Czoverall
            Ihor = (1/12)*tf2*math.pow(b2,3) + Ahor*math.pow(dz,2)

            Iz = Iver + Ihor

        return Iz

    def get_moment_of_interia_iacs(self, efficent_se=None, only_stf = False, tf1 = None):
        if only_stf:
            tf1 = 0
            b1 = 0
        else:
            tf1 = self._plate_th if tf1 == None else tf1
            b1 = self._spacing if efficent_se==None else efficent_se
        h = self._flange_th+self._web_height+tf1
        tw = self._web_th
        hw = self._web_height
        tf2 = self._flange_th
        b2 = self._flange_width

        Af = b2*tf2
        Aw = hw*tw

        ef = hw + tf2/2

        Iy = (Af*math.pow(ef,2)*math.pow(b2,2)/12) * ( (Af+2.6*Aw) / (Af+Aw))
        return Iy

    def get_torsional_moment_venant(self, reduced_tw = None, efficient_flange = True):
        # if efficient_flange:
        #     ef = self.get_ef_iacs()*1000
        # else:
        #     ef = self._flange_width * 1000
        tf = self._flange_th*1000
        tw = self._web_th*1000 if reduced_tw is None else reduced_tw
        bf = self._flange_width*1000
        hw = self._web_height*1000

        # if self._stiffener_type == 'FB':
        #     It = ((hw*math.pow(tw,3)/3e4) * (1-0.63*(tw/hw)) )
        # else:
        #     It = ((((ef-0.5*tf)*math.pow(tw,3))/3e4) * (1-0.63*(tw/(ef-0.5*tf))) + ((bf*math.pow(tf,3))/3e4)
        #           * (1-0.63*(tf/bf)) )
        # G = 80769.2
        # It2 = (2/3) * (math.pow(tw,3)*hw + bf*math.pow(tf, 3)) *(hw+tf/2)
        # print(It, It2*G)
        # print(hw, tw, bf, tf)
        I_t1 = 1.0 / 3.0 * math.pow(tw , 3) * hw + 1.0 / 3.0 * math.pow(tf, 3) * bf
        # I_t2 = 1.0 / 3.0 * math.pow(tw , 3) * (hw + tf) + 1.0 / 3.0 * math.pow(tf, 3) * (bf - tw)
        # print('It', I_t1, I_t2, It* 1e4)

        return I_t1#  * 1e4

    def get_polar_moment(self, reduced_tw  = None):
        tf = self._flange_th*1000
        tw = self._web_th*1000 if reduced_tw is None else reduced_tw
        ef = self.get_flange_eccentricity()*1000
        hw = self._web_height*1000
        b = self._flange_width*1000

        #Ipo = (A|w*(ef-0.5*tf)**2/3+Af*ef**2)*10e-4 #polar moment of interia in cm^4
        #Ipo = (tw/3)*math.pow(hw, 3) + tf*(math.pow(hw+tf/2,2)*b)+(tf/3)*(math.pow(ef+b/2,3)-math.pow(ef-b/2,3))

        # C24/3*C70^3+C26*((C70+C26/2)^2*C25)+C26/3*((C72+C25/2)^3-(C72-C25/2)^3) + (C25*C26^3)/12 + (C70*C24^3)/12
        Ipo = tw/3*math.pow(hw, 3)+tf*(math.pow(hw+tf/2,2)*b)+tf/3*(math.pow(ef+b/2,3)-math.pow(ef-b/2,3)) + \
              (b*math.pow(tf,3))/12 + (hw*math.pow(tw,3))/12

        return Ipo

    def get_flange_eccentricity(self):
        ecc = 0 if self._stiffener_type in ['FB', 'T'] else self._flange_width / 2 - self._web_th / 2
        return ecc

    def get_ef_iacs(self):

        if self._stiffener_type == 'FB':
            ef = self._web_height
        # elif self._stiffener_type == 'L-bulb':
        #     ef = self._web_height-0.5*self._flange_th
        elif self._stiffener_type in ['L', 'T', 'L-bulb', 'HP-profile', 'HP', 'HP-bulb', 'bulb']:
            ef = self._web_height + 0.5*self._flange_th
        return ef

    def get_stf_cog_eccentricity(self):
        e = (self._web_height * self._web_th * (self._web_height / 2) + self._flange_width * self._flange_th *
             (self._web_height + self._flange_th / 2)) / (self._web_height * self._web_th + self._flange_width * self._flange_th)
        return e

    def get_structure_prop(self):
        return self._main_dict

    def get_structure_type(self):
        return self._structure_type

    def get_stiffener_type(self):
        return self._stiffener_type

    def get_shear_area(self):
        '''
        Returning the shear area in [m^2]
        :return:
        '''
        return ((self._flange_th*self._web_th) + (self._web_th*self._plate_th) + (self._web_height*self._web_th))

    def set_main_properties(self, main_dict):
        '''
        Resettting all properties
        :param input_dictionary:
        :return:
        '''

        self._main_dict = main_dict
        self._plate_th = main_dict['plate_thk'][0]
        self._web_height = main_dict['stf_web_height'][0]
        self._web_th = main_dict['stf_web_thk'][0]
        self._flange_width = main_dict['stf_flange_width'][0]
        self._flange_th = main_dict['stf_flange_thk'][0]
        self._mat_yield = main_dict['mat_yield'][0]
        self._mat_factor = main_dict['mat_factor'][0]
        self._span = main_dict['span'][0]
        self._spacing = main_dict['spacing'][0]
        self._structure_type = main_dict['structure_type'][0]
        self._sigma_y1=main_dict['sigma_y1'][0]
        self._sigma_y2=main_dict['sigma_y2'][0]
        self._sigma_x1 = main_dict['sigma_x1'][0]
        self._sigma_x2 = main_dict['sigma_x2'][0]
        self._tauxy=main_dict['tau_xy'][0]
        self._plate_kpp = main_dict['plate_kpp'][0]
        self._stf_kps = main_dict['stf_kps'][0]
        self._km1 = main_dict['stf_km1'][0]
        self._km2 = main_dict['stf_km2'][0]
        self._km3 = main_dict['stf_km3'][0]
        self._stiffener_type=main_dict['stf_type'][0]
        try:
            self._girder_lg=main_dict['girder_lg'][0]
        except KeyError:
            self._girder_lg = 10
        try:
            self._pressure_side = main_dict['press_side'][0]
        except KeyError:
            self._pressure_side = 'p'
        self._zstar_optimization = main_dict['zstar_optimization'][0]
        self._puls_method = main_dict['puls buckling method'][0]
        self._puls_boundary = main_dict['puls boundary'][0]
        self._puls_stf_end  = main_dict['puls stiffener end'][0]
        self._puls_sp_or_up = main_dict['puls sp or up'][0]
        self._puls_up_boundary = main_dict['puls up boundary'][0]
        self._panel_or_shell = main_dict['panel or shell'][0]

    def set_stresses(self,sigy1,sigy2,sigx1,sigx2,tauxy):
        '''
        Setting the global stresses.
        :param sigy1:
        :param sigy2:
        :param sigx:
        :param tauxy:
        :return:
        '''
        self._main_dict['sigma_y1'][0]= sigy1
        self._sigma_y1 = sigy1

        self._main_dict['sigma_y2'][0]= sigy2
        self._sigma_y2  = sigy2

        self._main_dict['sigma_x1'][0]= sigx1
        self._sigma_x1 = sigx1

        self._main_dict['sigma_x2'][0]= sigx2
        self._sigma_x2 = sigx2

        self._main_dict['tau_xy'][0]= tauxy
        self._tauxy  = tauxy

    def get_cross_section_area(self, efficient_se = None, include_plate = True):
        '''
        Returns the cross section area.
        :return:
        '''
        tf1 = self._plate_th if include_plate else 0
        tf2 = self._flange_th
        if include_plate:
            b1 = self._spacing if efficient_se==None else efficient_se
        else:
            b1 = 0
        b2 = self._flange_width
        #h = self._flange_th+self._web_height+self._plate_th
        h = self._web_height
        tw = self._web_th
        #print('Plate: thk', tf1, 's', b1, 'Flange: thk', tf2, 'width', b2, 'Web: thk', tw, 'h', h)
        return tf1 * b1 + tf2 * b2 + h * tw

    def get_cross_section_centroid_with_effective_plate(self, se = None, tf1 = None, include_plate = True,
                                                        reduced_tw = None):
        '''
        Returns cross section centroid
        :return:
        '''
        # checked with example
        if include_plate:
            tf1 = self._plate_th if tf1 == None else tf1
            b1 = self._spacing if se == None else se
        else:
            tf1 = 0
            b1 = 0
        tf2 = self._flange_th
        b2 = self._flange_width
        tw = self._web_th if reduced_tw == None else reduced_tw/1000
        hw = self._web_height
        Ax = tf1 * b1 + tf2 * b2 + hw * tw
        effana = (tf1 * b1 * tf1/2 + hw * tw * (tf1 + hw / 2) + tf2 * b2 * (tf1+hw+tf2/2)) / Ax

        return effana

    def get_weight(self):
        '''
        Return the weight.
        :return:
        '''
        return 7850*self._span*(self._spacing*self._plate_th+self._web_height*self._web_th+self._flange_width*self._flange_th)

    def get_weight_width_lg(self):
        '''
        Return the weight including Lg
        :return:
        '''
        pl_area = self._girder_lg*self._plate_th
        stf_area = (self._web_height*self._web_th+self._flange_width*self._flange_th)*(self._girder_lg//self._spacing)
        return (pl_area+stf_area)*7850*self._span

    def set_span(self,span):
        '''
        Setting the span. Used when moving a point.
        :return: 
        '''
        self._span = span
        self._main_dict['span'][0] = span

    def get_puls_input(self, run_type: str = 'SP'):

        if self._stiffener_type == 'FB':
            stf_type = 'F'
        else:
            stf_type = self._stiffener_type
        map_boundary = {'Continuous': 'C', 'Sniped': 'S'}
        sig_x1 = self._sigma_x1
        sig_x2 = self._sigma_x2
        if sig_x1 * sig_x2 >= 0:
            sigxd = sig_x1 if abs(sig_x1) > abs(sig_x2) else sig_x2
        else:
            sigxd = max(sig_x1, sig_x2)
        puls_boundary = str(self._puls_boundary or 'Int').strip()
        axial_ml = 0 if puls_boundary in ('GT', 'Girder - trans') else sigxd
        trans_1_ml = 0 if puls_boundary in ('GL', 'Girder - long') else self._sigma_y1
        trans_2_ml = 0 if puls_boundary in ('GL', 'Girder - long') else self._sigma_y2
        if self._puls_sp_or_up == 'SP':
            return_dict = {'Identification': None, 'Length of panel': self._span*1000, 'Stiffener spacing': self._spacing*1000,
                            'Plate thickness': self._plate_th*1000,
                          'Number of primary stiffeners': 10,
                           'Stiffener type (L,T,F)': stf_type,
                            'Stiffener boundary': map_boundary[self._puls_stf_end]
                            if map_boundary[self._puls_stf_end] in ['C', 'S']
                            else 'C' if self._puls_stf_end == 'Continuous' else 'S',
                          'Stiff. Height': self._web_height*1000, 'Web thick.': self._web_th*1000,
                           'Flange width': self._flange_width*1000,
                            'Flange thick.': self._flange_th*1000, 'Tilt angle': 0,
                          'Number of sec. stiffeners': 0, 'Modulus of elasticity': 2.1e11/1e6, "Poisson's ratio": 0.3,
                          'Yield stress plate': self._mat_yield/1e6, 'Yield stress stiffener': self._mat_yield/1e6,
                            'Axial stress': 0 if self._puls_boundary == 'GT' else sigxd,
                           'Trans. stress 1': 0 if self._puls_boundary == 'GL' else self._sigma_y1,
                          'Trans. stress 2': 0 if self._puls_boundary == 'GL' else self._sigma_y2,
                           'Shear stress': self._tauxy,
                            'Pressure (fixed)': None, 'In-plane support': self._puls_boundary,
                           'sp or up': self._puls_sp_or_up}
        else:
            boundary = self._puls_up_boundary
            blist = list()
            if len(boundary) != 4:
                blist = ['SS', 'SS', 'SS', 'SS']
            else:
                for letter in boundary:
                    if letter.upper() == 'S':
                        blist.append('SS')
                    elif letter.upper() == 'C':
                        blist.append('CL')
                    else:
                        blist.append('SS')

            return_dict = {'Identification': None, 'Length of plate': self._span*1000, 'Width of c': self._spacing*1000,
                           'Plate thickness': self._plate_th*1000,
                         'Modulus of elasticity': 2.1e11/1e6, "Poisson's ratio": 0.3,
                          'Yield stress plate': self._mat_yield/1e6,
                         'Axial stress 1': 0 if self._puls_boundary == 'GT' else sigxd,
                           'Axial stress 2': 0 if self._puls_boundary == 'GT' else sigxd,
                           'Trans. stress 1': 0 if self._puls_boundary == 'GL' else self._sigma_y1,
                         'Trans. stress 2': 0 if self._puls_boundary == 'GL' else self._sigma_y2,
                           'Shear stress': self._tauxy, 'Pressure (fixed)': None, 'In-plane support': self._puls_boundary,
                         'Rot left': blist[0], 'Rot right': blist[1], 'Rot upper': blist[2], 'Rot lower': blist[3],
                           'sp or up': self._puls_sp_or_up}
        return return_dict

    def get_buckling_ml_input(self, design_lat_press: float = 0, sp_or_up: str = 'SP', alone = True, csr = False):
        '''
        Classes in data from ML

        {'negative utilisation': 1, 'non-zero': 2, 'Division by zero': 3, 'Overflow': 4, 'aspect ratio': 5,
        'global slenderness': 6, 'pressure': 7, 'web-flange-ratio': 8,  'below 0.87': 9,
                  'between 0.87 and 1': 10, 'above 1': 11}
        '''
        stf_type = {'T-bar': 1,'T': 1,  'L-bulb': 2, 'Angle': 3, 'Flatbar': 4, 'FB': 4, 'L': 3}
        stf_end = {'Cont': 1, 'C':1 , 'Sniped': 2, 'S': 2}
        field_type = {'Integrated': 1,'Int': 1, 'Girder - long': 2,'GL': 2, 'Girder - trans': 3,  'GT': 3}
        up_boundary = {'SS': 1, 'CL': 2}
        map_boundary = {'Continuous': 'C', 'Sniped': 'S'}
        sig_x1 = self._sigma_x1
        sig_x2 = self._sigma_x2
        if sig_x1 * sig_x2 >= 0:
            sigxd = sig_x1 if abs(sig_x1) > abs(sig_x2) else sig_x2
        else:
            sigxd = max(sig_x1, sig_x2)
        puls_boundary = str(self._puls_boundary or 'Int').strip()
        axial_ml = 0 if puls_boundary in ('GT', 'Girder - trans') else sigxd
        trans_1_ml = 0 if puls_boundary in ('GL', 'Girder - long') else self._sigma_y1
        trans_2_ml = 0 if puls_boundary in ('GL', 'Girder - long') else self._sigma_y2
        if self._puls_sp_or_up == 'SP':

            if csr == False:

                this_field =  [self._span * 1000, self._spacing * 1000, self._plate_th * 1000, self._web_height * 1000,
                               self._web_th * 1000, self._flange_width * 1000, self._flange_th * 1000, self._mat_yield / 1e6,
                               self._mat_yield / 1e6, axial_ml, trans_1_ml, trans_2_ml, self._tauxy,
                               design_lat_press/1000, stf_type[self._stiffener_type],
                               stf_end[map_boundary[self._puls_stf_end]]]
            else:
                this_field =  [self._span * 1000, self._spacing * 1000, self._plate_th * 1000, self._web_height * 1000,
                               self._web_th * 1000, self._flange_width * 1000, self._flange_th * 1000, self._mat_yield / 1e6,
                               self._mat_yield / 1e6,  axial_ml, trans_1_ml, trans_2_ml, self._tauxy,
                               design_lat_press/1000, stf_type[self._stiffener_type],
                               stf_end[map_boundary[self._puls_stf_end]],
                               field_type[puls_boundary]]
        else:
            boundary = str(self._puls_up_boundary or 'SSSS').strip().upper()
            if len(boundary) != 4:
                boundary = 'SSSS'
            ss_cl_list = list()
            for letter_i in boundary:
                if letter_i == 'C':
                    ss_cl_list.append(up_boundary['CL'])
                else:
                    ss_cl_list.append(up_boundary['SS'])
            b1, b2, b3, b4 = ss_cl_list
            if csr == False:
                this_field =  [self._span * 1000, self._spacing * 1000, self._plate_th * 1000, self._mat_yield / 1e6,
                               axial_ml, trans_1_ml, trans_2_ml, self._tauxy, design_lat_press/1000,
                               b1, b2, b3, b4]
            else:
                this_field =  [self._span * 1000, self._spacing * 1000, self._plate_th * 1000, self._mat_yield / 1e6,
                               axial_ml, trans_1_ml, trans_2_ml, self._tauxy, design_lat_press/1000,
                               field_type[puls_boundary], b1, b2, b3, b4]
        if alone:
            return [this_field,]
        else:
            return this_field

class CalcScantlings(Structure):
    '''
    This Class does the calculations for the plate fields. 
    Input is a structure object, same as for the structure class.
    The class inherits from Structure class.
    '''

    def __init__(self, main_dict: dict = None, lat_press = True, category = 'secondary'):
        super(CalcScantlings,self).__init__(main_dict=main_dict)

        self.lat_press = lat_press
        self.category = category
        self._need_recalc = True

    @property
    def need_recalc(self):
        return self._need_recalc

    @need_recalc.setter
    def need_recalc(self, val):
        self._need_recalc = val

    def get_results_for_report(self,lat_press=0):
        '''
        Returns a string for the report.
        :return:
        '''
        buckling = AllStructure(Plate=self, Stiffener=self, calculation_domain='Flat plate, stiffened')
        buckling.mat_yield = self._mat_yield
        buckling.lat_press = lat_press / 1000
        buc = buckling.plate_buckling()
        min_section_modulus = self.get_dnv_min_section_modulus(design_pressure_kpa=lat_press)
        min_section_modulus_text = 'inf' if math.isinf(min_section_modulus) else \
            str(int(min_section_modulus*1000**3))

        return 'Minimum section modulus:'\
               +min_section_modulus_text\
               +'mm^3 '+' Minium plate thickness: '\
               +str(round(self.get_dnv_min_thickness(design_pressure_kpa=lat_press),1))+\
               ' Buckling results: plate: '+str(round(buc['Plate']['Plate buckling'], 1))+\
               ' stiffener plate side: '+str(round(buc['Stiffener']['Overpressure plate side'], 1))+\
               ' stiffener side: '+str(round(buc['Stiffener']['Overpressure stiffener side'], 1))

    def calculate_slamming_plate(self, slamming_pressure, red_fac = 1):
        ''' Slamming pressure input is Pa '''
        ka1 = 1.1
        ka2 = min(max(0.4, self._spacing / self._span), 1)

        ka = math.pow(ka1 - 0.25*ka2,2)
        sigmaf = self._mat_yield/1e6  # MPa

        psl = red_fac * slamming_pressure/1000  # kPa
        Cd = 1.5

        return 0.0158*ka*self._spacing*1000*math.sqrt(psl/(Cd*sigmaf))

    def calculate_slamming_stiffener(self, slamming_pressure, angle = 90, red_fac = 1):
        tk = 0
        psl = slamming_pressure / 1000  # kPa
        Pst = psl * red_fac  # Currently DNV does not use psl/2 for slamming.
        sigmaf = self._mat_yield / 1e6  # MPa
        hw, twa, tp, tf, bf, s = [(val - tk) * 1000 for val in [self._web_height, self._web_th, self._plate_th,
                                                                self._flange_th, self._flange_width, self._spacing]]
        ns = 2
        tau_eH = sigmaf/math.sqrt(3)
        h_stf = (self._web_height+self._flange_th)*1000
        f_shr = 0.7
        lbdg = self._span
        lshr = self._span - self._spacing/4000
        dshr = h_stf + tp if 75 <= angle <= 90 else (h_stf + tp)*math.sin(math.radians(angle))
        tw = (f_shr*Pst*s*lshr)/(dshr*tau_eH)

        if self._web_th*1000 < tw:
            return {'tw_req': tw, 'Zp_req':None}
        fpl = 8* (1+(ns/2))
        Zp_req = (1.2*Pst*s*math.pow(lbdg,2)/(fpl*sigmaf)) + \
                  (ns*(1-math.sqrt(1-math.pow(tw/twa,2)))*hw*tw*(hw+tp))/8000

        return {'tw_req': tw, 'Zp_req':Zp_req}

    def check_all_slamming(self, slamming_pressure, stf_red_fact = 1, pl_red_fact = 1):
        ''' A summary check of slamming '''

        pl_chk = self.calculate_slamming_plate(slamming_pressure, red_fac= pl_red_fact)
        if self._plate_th*1000 < pl_chk:
            chk1 = pl_chk / (self._plate_th*1000)
            return False, chk1

        stf_res = self.calculate_slamming_stiffener(slamming_pressure, red_fac = stf_red_fact)
        #print('Slamming checked')
        if self._web_th*1000 < stf_res['tw_req']:
            chk2 = stf_res['tw_req'] / (self._web_th*1000)
            return False, chk2

        if stf_res['Zp_req'] is not None:
            eff_pl_sec_mod = self.get_net_effective_plastic_section_modulus()
            if eff_pl_sec_mod < stf_res['Zp_req']:
                chk3 = stf_res['Zp_req']/eff_pl_sec_mod
                return False, chk3

        return True, None

    def get_net_effective_plastic_section_modulus(self, angle = 90):
        ''' Calculated according to Rules for classification: Ships — DNVGL-RU-SHIP Pt.3 Ch.3. Edition July 2017,
            page 83 '''
        tk = 0
        angle_rad = math.radians(angle)
        hw, tw, tp, tf, bf = [(val - tk) * 1000 for val in [self._web_height, self._web_th, self._plate_th, self._flange_th,
                                                            self._flange_width]]
        h_stf = (self._web_height+self._flange_th)*1000
        de_gr = 0
        tw_gr = self._web_th*1000
        hf_ctr = h_stf-0.5*tf if self.get_stiffener_type() not in ['L','L-bulb'] else h_stf - de_gr - 0.5*tf
        bf_ctr = 0 if self.get_stiffener_type() == 'T' else 0.5*(tf - tw_gr)
        beta = 0.5
        gamma = (1 + math.sqrt(3+12*beta))/4

        Af = 0 if self.get_stiffener_type() == 'FB' else bf*tf

        if 75 <= angle <= 90:
            zpl = (hw*tw*(hw+tp)/2000) + ( (2*gamma-1) * Af * ((hf_ctr + tp/2)) / 1000)
        elif angle < 75:
            zpl = (hw*tw*(hw+tp)/2000)+\
                  ( (2*gamma-1) * Af * ((hf_ctr + tp/2) * math.sin(angle_rad) - bf_ctr*math.cos(angle_rad)) / 1000)

        return zpl

    def get_dnv_min_section_modulus(self, design_pressure_kpa, printit = False):
        ''' Section modulus according to DNV rules '''

        design_pressure = design_pressure_kpa
        fy = self._mat_yield / 1e6
        fyd = fy/self._mat_factor

        sigma_y = self._sigma_y2 + (self._sigma_y1-self._sigma_y2)\
                                       *(min(0.25*self._span,0.5*self._spacing)/self._span)
        sig_x1 = self._sigma_x1
        sig_x2 = self._sigma_x2
        if sig_x1 * sig_x2 >= 0:
            sigxd = sig_x1 if abs(sig_x1) > abs(sig_x2) else sig_x2
        else:
            sigxd =max(sig_x1 , sig_x2)

        sigma_jd = math.sqrt(math.pow(sigxd,2)+math.pow(sigma_y,2)-
                             sigxd*sigma_y+3*math.pow(self._tauxy,2))

        sigma_pd2 = fyd-sigma_jd  # design_bending_stress_mpa
        if sigma_pd2 <= 0:
            return math.inf

        kps = self._stf_kps  # 1 is clamped, 0.9 is simply supported.
        km_sides = min(self._km1,self._km3)  # see table 3 in DNVGL-OS-C101 (page 62)
        km_middle = self._km2  # see table 3 in DNVGL-OS-C101 (page 62)

        Zs = ((math.pow(self._span, 2) * self._spacing * design_pressure) /
              (min(km_middle, km_sides) * (sigma_pd2) * kps)) * math.pow(10, 6)
        if printit:
            print('Sigma y1', self._sigma_y1, 'Sigma y2', self._sigma_y2, 'Sigma x', self._sigma_x1,
                  'Pressure', design_pressure, 'fy', fy,
                  'Section mod', max(math.pow(15, 3) / math.pow(1000, 3), Zs / math.pow(1000, 3)))
        return max(math.pow(15, 3) / math.pow(1000, 3), Zs / math.pow(1000, 3))

    def get_dnv_min_thickness(self, design_pressure_kpa):
        '''
        Return minimum thickness in mm
        :param design_pressure_kpa:
        :return:
        '''

        design_pressure = design_pressure_kpa
        #print(self._sigma_x1)
        self.span
        sigma_y = self._sigma_y2 + (self._sigma_y1-self._sigma_y2)\
                                       *(min(0.25*self.span,0.5*self._spacing)/self._span)

        sig_x1 = self._sigma_x1
        sig_x2 = self._sigma_x2
        if sig_x1 * sig_x2 >= 0:
            sigxd = sig_x1 if abs(sig_x1) > abs(sig_x2) else sig_x2
        else:
            sigxd =max(sig_x1 , sig_x2)

        sigma_jd = math.sqrt(math.pow(sigxd,2)+math.pow(sigma_y,2)-
                             sigxd*sigma_y+3*math.pow(self._tauxy,2))

        fy = self._mat_yield / 1000000
        fyd = fy/self._mat_factor
        if fyd - sigma_jd <= 0:
            return math.inf
        sigma_pd1 = min(1.3*(fyd-sigma_jd), fyd)
        #print(fyd, sigma_jd, fyd)
        if self.category == 'secondary':
            t0 = 5
        else:
            t0 = 7

        t_min = (14.3 * t0) / math.sqrt(fyd)

        ka = math.pow(1.1 - 0.25  * self._spacing/self._span, 2)

        if ka > 1:
            ka =1
        elif ka<0.72:
            ka = 0.72

        assert sigma_pd1 > 0, 'sigma_pd1 must be positive | current value is: ' + str(sigma_pd1)
        assert self._plate_kpp is not None, 'Fixation parameters must be set.'
        t_min_bend = (15.8 * ka * self._spacing * math.sqrt(design_pressure)) / \
                     math.sqrt(sigma_pd1 *self._plate_kpp)

        if self.lat_press:
            return max(t_min, t_min_bend)
        else:
            return t_min

    def get_minimum_shear_area(self, pressure):
        '''
        Calculating minimum section area according to ch 6.4.4.

        Return [m^2]
        :return:
        '''
        #print('SIGMA_X ', self._sigma_x1)
        l = self._span
        s = self._spacing
        fy = self._mat_yield

        fyd = (fy/self._mat_factor)/1e6 #yield strength
        sig_x1 = self._sigma_x1
        sig_x2 = self._sigma_x2
        if sig_x1 * sig_x2 >= 0:
            sigxd = sig_x1 if abs(sig_x1) > abs(sig_x2) else sig_x2
        else:
            sigxd =max(sig_x1 , sig_x2)

        shear_reserve = math.pow(fyd, 2) - math.pow(sigxd, 2)
        if shear_reserve <= 0:
            return math.inf
        taupds = 0.577*math.sqrt(shear_reserve)

        As = ((l*s*pressure)/(2*taupds)) * math.pow(10,3)

        return As/math.pow(1000,2)

    def is_acceptable_sec_mod(self, section_module, pressure):
        '''
        Checking if the result is accepable.
        :param section_module:
        :param pressure:
        :return:
        '''

        return min(section_module) >= self.get_dnv_min_section_modulus(pressure)

    def is_acceptable_shear_area(self, shear_area, pressure):
        '''
        Returning if the shear area is ok.
        :param shear_area:
        :param pressure:
        :return:
        '''

        return shear_area >= self.get_minimum_shear_area(pressure)

    def get_plate_efficent_b(self,design_lat_press=0,axial_stress=50,
                                 trans_stress_small=100,trans_stress_large=100):
        '''
        Simple buckling calculations according to DNV-RP-C201
        :return:
        '''

        #7.2 Forces in the idealised stiffened plate

        s = self._spacing #ok
        t = self._plate_th #ok
        l = self._span #ok

        E = 2.1e11 #ok

        pSd = design_lat_press*1000
        sigy1Sd =trans_stress_large*1e6
        sigy2Sd =trans_stress_small*1e6
        sigxSd = axial_stress*1e6

        fy = self._mat_yield #ok

        #7.3 Effective plate width
        alphap=0.525*(s/t)*math.sqrt(fy/E) # reduced plate slenderness, checked not calculated with ex
        alphac = 1.1*(s/t)*math.sqrt(fy/E) # checked not calculated with example
        mu6_9 = 0.21*(alphac-0.2)

        if alphac<=0.2: kappa = 1 # eq6.7, all kappa checked not calculated with example
        elif 0.2<alphac<2: kappa = (1/(2*math.pow(alphac,2)))*(1+mu6_9+math.pow(alphac,2)
                                                               -math.sqrt(math.pow(1+mu6_9+math.pow(alphac,2),2)
                                                                          -4*math.pow(alphac,2))) # ok
        else: kappa=(1/(2*math.pow(alphac,2)))+0.07 # ok

        ha = max(0, 0.05*(s/t)-0.75)
        kp = 1 if pSd<=2*((t/s)**2)*fy else 1-ha*((pSd/fy)-2*(t/s)**2)

        sigyR=( (1.3*t/l)*math.sqrt(E/fy)+kappa*(1-(1.3*t/l)*math.sqrt(E/fy)))*fy*kp # checked not calculated with example
        l1 = min(0.25*l,0.5*s)

        sig_min, sig_max = min(sigy1Sd,sigy2Sd),max(sigy1Sd,sigy2Sd) # self-made
        sigySd = sig_min+(sig_max-sig_min)*(1-l1/l) # see 6.8, page 15

        ci = 1-s/(120*t) if (s/t)<=120 else 0 # checked not calculated with example

        Cxs = (alphap-0.22)/math.pow(alphap,2) if alphap > 0.673 else 1 # reduction factor longitudinal
        # eq7.16, reduction factor transverse, compression (positive) else tension

        if sigySd >= 0:
            cys_term = 1-math.pow(sigySd/sigyR,2) + ci*((sigxSd*sigySd)/(Cxs*fy*sigyR)) if Cxs*fy*sigyR != 0 else 0
            Cys = 0 if cys_term <= 0 else math.sqrt(cys_term)
        else:
            cys_term = 4-3*math.pow(sigySd/fy,2) if fy != 0 else 0
            Cys = min(0.5*((0 if cys_term <= 0 else math.sqrt(cys_term))+sigySd/fy),1) #ok, checked
            Cys = max(0, Cys)

        #7.7.3 Resistance parameters for stiffeners
        return s * Cxs * Cys # 7.3, eq7.13, che

    def buckling_local_stiffener(self):
        '''
        Local requirements for stiffeners. Chapter 9.11.
        :return:
        '''

        epsilon = math.sqrt(235 / (self._mat_yield / 1e6))

        if self._stiffener_type in ['L', 'L-bulb']:
            c = self._flange_width - self._web_th/2
        elif self._stiffener_type == 'T':
            c = self._flange_width/2 - self._web_th/2
        elif self._stiffener_type == 'FB':
            return self._web_height <= 42 * self._web_th * epsilon, self._web_height/(42 * self._web_th * epsilon)

        # print(self._web_height, self._web_th, self._flange_width ,self._flange_th )
        # print('c:',c, 14 * self._flange_th * epsilon, ' | ',  self._web_height, 42 * self._web_th * epsilon)
        # print(c <= (14  * self._flange_th * epsilon) and self._web_height <= 42 * self._web_th * epsilon)
        # print(c/(14  * self._flange_th * epsilon), self._web_height / (42 * self._web_th * epsilon))
        # print('')

        return c <= (14  * self._flange_th * epsilon) and self._web_height <= 42 * self._web_th * epsilon, \
               max(c/(14  * self._flange_th * epsilon), self._web_height / (42 * self._web_th * epsilon))

    def is_acceptable_pl_thk(self, design_pressure):
        '''
        Checking if the thickness is acceptable.
        :return:
        '''
        return self.get_dnv_min_thickness(design_pressure) <= self._plate_th*1000

class AllStructure():
    '''
    Calculation of structure
    '''
    def __init__(self, Plate: CalcScantlings = None, Stiffener: CalcScantlings = None, Girder: CalcScantlings = None,
                 main_dict = None, calculation_domain: str = None):
        super(AllStructure, self).__init__()
        self._Plate = Plate  # This contain the stresses
        self._Stiffener = Stiffener
        self._Girder = Girder
        self._lat_press = 0

        self._v = 0.3
        self._E = 2.1e11
        if main_dict is not None:
            self._min_lat_press_adj_span = None if main_dict['minimum pressure in adjacent spans'][0] == 0 else \
                main_dict['minimum pressure in adjacent spans'][0]
            self._mat_yield =  main_dict['material yield'][0]
            self._stress_load_factor = main_dict['load factor on stresses'][0]
            self._lat_load_factor = main_dict['load factor on pressure'][0]
            self._method = main_dict['buckling method'][0]
            self._stf_end_support = main_dict['stiffener end support'][0]#'Continuous'
            self._girder_end_support = main_dict['girder end support'][0]#'Continuous'
            self._tension_field_action = main_dict['tension field'][0]# 'not allowed'
            self._stiffened_plate_effective_aginst_sigy = main_dict['plate effective agains sigy'][0] #True
            self._buckling_length_factor_stf = None if main_dict['buckling length factor stf'][0] == 0 else \
                main_dict['buckling length factor stf'][0]
            self._buckling_length_factor_girder = None if main_dict['buckling length factor girder'][0] == 0 else \
                main_dict['buckling length factor girder'][0]
            self._km3 = main_dict['km3'][0]#12
            self._km2 = main_dict['km2'][0]#24
            self._stf_dist_between_lateral_supp = None if main_dict['stiffener distance between lateral support'][0] == 0 \
                else main_dict['stiffener distance between lateral support'][0]
            self._girder_dist_between_lateral_supp = None if main_dict['girder distance between lateral support'][0] == 0 \
                else main_dict['girder distance between lateral support'][0]
            self._panel_length_Lp = None if main_dict['panel length, Lp'][0] == 0 else main_dict['panel length, Lp'][0]
            self._overpressure_side = main_dict['pressure side'][0] # either 'stiffener side', 'plate side', 'both sides'
            self._fab_method_stiffener = main_dict['fabrication method stiffener'][0]#'welded'
            self._fab_method_girder = main_dict['fabrication method girder'][0]#'welded'
            self._calculation_domain = main_dict['calculation domain'][0]
            self._need_recalc = True
        else: # setting default values
            self._stf_end_support = 'Continuous'
            self._girder_end_support = 'Continuous'
            self._tension_field_action = 'not allowed'
            self._stiffened_plate_effective_aginst_sigy = ['Stf. pl. effective against sigma y',
                                                          'All sigma y to girder'][0]
            self._km3 = 12
            self._km2 = 24
            self._overpressure_side = 'both sides'
            self._fab_method_stiffener = 'welded'
            self._fab_method_girder = 'welded'
            self._need_recalc = True
            self._mat_yield = None
            self._lat_load_factor = 1
            self._stress_load_factor = 1
            self._lat_press = 0
            self._min_lat_press_adj_span = None
            self._buckling_length_factor_stf = None
            self._buckling_length_factor_girder = None
            self._stf_dist_between_lateral_supp = None
            self._girder_dist_between_lateral_supp = None
            self._panel_length_Lp = None
            self._calculation_domain = calculation_domain
            
    
    @property
    def method(self):
        return self._method
    @method.setter
    def method(self, val):
        self._method = val
    @property
    def tension_field_action(self):
        return self._tension_field_action
    @tension_field_action.setter
    def tension_field_action(self, val):
        self._tension_field_action = val
    @property
    def stiffened_plate_effective_aginst_sigy(self):
        return self._stiffened_plate_effective_aginst_sigy
    @stiffened_plate_effective_aginst_sigy.setter
    def stiffened_plate_effective_aginst_sigy(self, val):
        self._stiffened_plate_effective_aginst_sigy = val
    @property
    def km3(self):
        return self._km3
    @km3.setter
    def km3(self, val):
        self._km3 = val
    @property
    def km2(self):
        return self._km2
    @km2.setter
    def km2(self, val):
        self._km2 = val
    
    @property
    def need_recalc(self):
        return self._need_recalc
    @need_recalc.setter
    def need_recalc(self, val):
        self._need_recalc = val
    
    @property
    def overpressure_side(self):
        return self._overpressure_side
    @overpressure_side.setter
    def overpressure_side(self, val):
        self._overpressure_side = val
    @property
    def fab_method_stiffener(self):
        return self._fab_method_stiffener
    @fab_method_stiffener.setter
    def fab_method_stiffener(self, val):
        self._fab_method_stiffener = val
    
    @property
    def fab_method_girder(self):
        return self._fab_method_girder
    @fab_method_girder.setter
    def fab_method_girder(self, val):
        self._fab_method_girder = val
    @property
    def stf_end_support(self):
        return self._stf_end_support
    @stf_end_support.setter
    def stf_end_support(self, val):
        self._stf_end_support = val
    @property
    def girder_end_support(self):
        return self._girder_end_support
    @girder_end_support.setter
    def girder_end_support(self, val):
        self._girder_end_support = val
    @property
    def need_recalc(self):
        return self._need_recalc
    @need_recalc.setter
    def need_recalc(self, val):
        self._need_recalc = val
        
    @property
    def lat_press(self):
        return self._lat_press

    @lat_press.setter
    def lat_press(self, val):
        self._lat_press = val
        
    @property
    def Plate(self):
        return self._Plate

    @Plate.setter
    def Plate(self, val):
        self._Plate = val
        
    @property
    def Stiffener(self):
        return self._Stiffener
    @Stiffener.setter
    def Stiffener(self, val):
        self._Stiffener = val
        
    @property
    def Girder(self):
        return self._Girder
    @Girder.setter
    def Girder(self, val):
        self._Girder = val
    
    @property
    def overpressure_side(self):
        return self._overpressure_side
    @overpressure_side.setter
    def overpressure_side(self, val):
        self._overpressure_side = val
    
    @property
    def calculation_domain(self):
        return self._calculation_domain
    @calculation_domain.setter
    def calculation_domain(self, val):
        self._calculation_domain = val

    @property
    def E(self):
        assert self._E is not None, 'Missing variable: self._E'
        return self._E
    @E.setter
    def E(self, val):
        self._E = val
    @property
    def v(self):
        assert self._v is not None, 'Missing variable: self._v'
        return self._v
    @v.setter
    def v(self, val):
        self._v = val
    @property
    def mat_yield(self):
        assert self._mat_yield is not None, 'Missing variable: self._mat_yield'
        return self._mat_yield
    @mat_yield.setter
    def mat_yield(self, val):
        self._mat_yield = val

    def get_method(self):
        gird_opt = ['Stf. pl. effective against sigma y', 'All sigma y to girder']
        #stf_opt = ['allowed', 'not allowed']
        # if self.calculation_domain == "Flat plate, stiffened with girder":

        if self._stiffened_plate_effective_aginst_sigy == True:
            self._stiffened_plate_effective_aginst_sigy = gird_opt[0]
        elif self._stiffened_plate_effective_aginst_sigy == False:
            self._stiffened_plate_effective_aginst_sigy = gird_opt[1]

        if self.calculation_domain == "Flat plate, stiffened with girder":
            if self._stiffened_plate_effective_aginst_sigy == gird_opt[0]:
                return 1
            else:
                return 2
        else:
            return 1
        # else:
        #     if self._tension_field_action == stf_opt[0]:
        #         return 1
        #     else:
        #         return 2

    def get_main_properties(self):
        main_dict = dict()
        main_dict['minimum pressure in adjacent spans'] = [self._min_lat_press_adj_span,  '']
        main_dict['material yield'] = [self._mat_yield, 'Pa']
        main_dict['load factor on stresses'] = [self._stress_load_factor, '']
        main_dict['load factor on pressure'] = [self._lat_load_factor, '']
        main_dict['buckling method'] = [self._method, '']
        main_dict['stiffener end support'] = [self._stf_end_support, '']  # 'Continuous'
        main_dict['girder end support'] = [self._girder_end_support, '']  # 'Continuous'
        main_dict['tension field'] = [self._tension_field_action, '']  # 'not allowed'
        main_dict['plate effective agains sigy'] = [self._stiffened_plate_effective_aginst_sigy, '']  # True
        main_dict['buckling length factor stf'] = [self._buckling_length_factor_stf, '']
        main_dict['buckling length factor girder'] = [self._buckling_length_factor_girder, '']
        main_dict['km3'] = [self._km3, '']  # 12
        main_dict['km2'] = [self._km2, '']  # 24
        main_dict['girder distance between lateral support'] = [self._girder_dist_between_lateral_supp, '']
        main_dict['stiffener distance between lateral support'] = [self._stf_dist_between_lateral_supp, '']
        main_dict['panel length, Lp'] = [self._panel_length_Lp, '']
        main_dict['pressure side'] = [self._overpressure_side, '']  # either 'stiffener', 'plate', 'both'
        main_dict['fabrication method stiffener'] = [self._fab_method_stiffener, '']
        main_dict['fabrication method girder'] = [self._fab_method_girder, '']
        main_dict['calculation domain']= [self._calculation_domain, '']

        return {'main dict': main_dict, 'Plate': self._Plate.get_structure_prop(),
                'Stiffener': None if self._Stiffener is None else self._Stiffener.get_structure_prop(),
                'Girder': None if self._Girder is None else self._Girder.get_structure_prop()}

    def set_main_properties(self, prop_dict):
        main_dict = prop_dict['main dict']
        self._min_lat_press_adj_span = None if main_dict['minimum pressure in adjacent spans'][0] == 0 else \
            main_dict['minimum pressure in adjacent spans'][0]
        self._mat_yield =  main_dict['material yield'][0]
        self._stress_load_factor = main_dict['load factor on stresses'][0]
        self._lat_load_factor = main_dict['load factor on pressure'][0]
        self._method = main_dict['buckling method'][0]
        self._stf_end_support = main_dict['stiffener end support'][0]#'Continuous'
        self._girder_end_support = main_dict['girder end support'][0]#'Continuous'
        self._tension_field_action = main_dict['tension field'][0]# 'not allowed'
        self._stiffened_plate_effective_aginst_sigy = main_dict['plate effective agains sigy'][0] #True
        self._buckling_length_factor_stf = None if main_dict['buckling length factor stf'][0] == 0 else \
            main_dict['buckling length factor stf'][0]
        self._buckling_length_factor_girder = None if main_dict['buckling length factor girder'][0] == 0 else \
            main_dict['buckling length factor girder'][0]
        self._km3 = main_dict['km3'][0]#12
        self._km2 = main_dict['km2'][0]#24
        self._girder_dist_between_lateral_supp = None if main_dict['girder distance between lateral support'][0] in [0, None, ''] else \
            main_dict['girder distance between lateral support'][0]
        self._stf_dist_between_lateral_supp = None if main_dict['stiffener distance between lateral support'][0]  in [0, None, ''] else \
            main_dict['stiffener distance between lateral support'][0]
        self._panel_length_Lp = None if main_dict['panel length, Lp'][0] == 0 else main_dict['panel length, Lp'][0]
        self._overpressure_side = main_dict['pressure side'][0] # either 'stiffener', 'plate', 'both'
        self._fab_method_stiffener = main_dict['fabrication method stiffener'][0]#'welded'
        self._fab_method_girder = main_dict['fabrication method girder'][0]#'welded'

        self._Plate.set_main_properties(prop_dict['Plate'])

        if prop_dict['Stiffener'] is not None and self._Stiffener is None:
            self._Stiffener = CalcScantlings(prop_dict['Stiffener'])
        elif prop_dict['Stiffener'] is not None and self._Stiffener is not None:
            self._Stiffener.set_main_properties(prop_dict['Stiffener'])
        else:
            self._Stiffener = None

        if prop_dict['Girder'] is not None and self._Girder is None:
            self._Girder = CalcScantlings(prop_dict['Girder'])
        elif prop_dict['Girder'] is not None and self._Girder is not None:
            self._Girder.set_main_properties(prop_dict['Girder'])
        else:
            self._Girder = None

        self._calculation_domain = main_dict['calculation domain'][0]

    def plate_buckling(self, optimizing = False):
        '''
        Summary
        '''
        return_dummy = {'Plate': {'Plate buckling': 0},
                        'Stiffener': {'Overpressure plate side': 0, 'Overpressure stiffener side': 0,
                                      'Resistance between stiffeners': 0, 'Shear capacity': 0},
                        'Girder': {'Overpressure plate side': 0, 'Overpressure girder side': 0, 'Shear capacity': 0},
                        'Local buckling': 0}

        unstf_pl = self.unstiffened_plate_buckling(optimizing = optimizing)
        up_buckling = max([unstf_pl['UF Pnt. 5  Lateral loaded plates'], unstf_pl['UF sjsd'],
                           max([unstf_pl['UF Longitudinal stress'],  unstf_pl['UF transverse stresses'],
                                unstf_pl['UF Shear stresses'], unstf_pl['UF Combined stresses']])
                           if all([self._Girder is None, self._Stiffener is None]) else 0])
        if optimizing and up_buckling > 1:
            return_dummy['Plate']['Plate buckling'] = up_buckling
            return return_dummy

        if not optimizing:
            local_buckling = self.local_buckling()

        if self._Stiffener is not None:
            stf_pla = self.stiffened_panel(unstf_pl_data=unstf_pl, optimizing=optimizing)
            if all([optimizing, type(stf_pla) == list]):
                return_dummy['Stiffener'][stf_pla[0]] = stf_pla[1]
                return return_dummy

            stf_buckling_pl_side = stf_pla['UF Plate side'] if self._stf_end_support == 'Continuous' else \
                stf_pla['UF simply supported plate side']
            stf_buckling_stf_side = stf_pla['UF Stiffener side'] if self._stf_end_support == 'Continuous' else \
                stf_pla['UF simply supported stf side']
            stf_plate_resistance = stf_pla['UF Plate resistance']
            stf_shear_capacity = stf_pla['UF Shear force']
        else:
            stf_buckling_pl_side, stf_buckling_pl_side, stf_buckling_stf_side, stf_plate_resistance, \
            stf_shear_capacity = 0,0,0,0,0

        if self._Girder is not None:
            girder = self.girder(unstf_pl_data=unstf_pl, stf_pl_data=stf_pla, optmizing=optimizing)
            if all([optimizing, type(girder) == list]):
                return_dummy['Girder'][girder[0]] = girder[1]
                return return_dummy

            girder_buckling_pl_side = girder['UF Cont. plate side'] if self._girder_end_support == 'Continuous' else \
                girder['UF Simplified plate side']
            girder_buckling_girder_side = girder['UF Cont. girder side'] if self._girder_end_support == 'Continuous' \
                else \
                girder['UF Simplified girder side']
            girder_shear_capacity = girder['UF shear force']
        else:
            girder_buckling_pl_side, girder_buckling_girder_side, girder_shear_capacity = 0,0,0
        
        return {'Plate': {'Plate buckling': up_buckling},
                'Stiffener': {'Overpressure plate side': stf_buckling_pl_side,
                                                    'Overpressure stiffener side': stf_buckling_stf_side, 
                                                    'Resistance between stiffeners': stf_plate_resistance,
                                                    'Shear capacity': stf_shear_capacity},
                'Girder': {'Overpressure plate side': girder_buckling_pl_side,
                           'Overpressure girder side': girder_buckling_girder_side,
                           'Shear capacity': girder_shear_capacity},
                'Local buckling': 0 if optimizing else local_buckling}

    def unstiffened_plate_buckling(self, optimizing = False):

        unstf_pl_data = dict()

        E = self._E/1e6
        v = self._v
        fy = self._mat_yield/1e6
        gammaM = self._Plate.mat_factor
        t = self._Plate.t
        s = self._Plate.spacing
        l = self._Plate.span*1000

        tsd = abs(self._Plate.tau_xy * self._stress_load_factor)
        psd = self._lat_press*self._lat_load_factor

        sig_x1 = self._Plate.sigma_x1*self._stress_load_factor
        sig_x2 = self._Plate.sigma_x2*self._stress_load_factor

        sig_y1 = self._Plate.sigma_y1 * self._stress_load_factor
        sig_y2 = self._Plate.sigma_y2 * self._stress_load_factor


        if sig_x1 * sig_x2 >= 0:
            Use_Smax_x = sxsd = sig_x1 if abs(sig_x1) > abs(sig_x2) else sig_x2
        else:
            Use_Smax_x = sxsd =max(sig_x1 , sig_x2)

        if sig_y1 * sig_y2 >= 0:
            Use_Smax_y = sy1sd = sig_y1 if abs(sig_y1) > abs(sig_y2) else sig_y2
        else:
            Use_Smax_y = sy1sd = max(sig_y1 , sig_y2)

        if sig_x1 * sig_x2 >= 0:
            Use_Smin_x = sig_x2 if abs(sig_x1) > abs(sig_x2) else sig_x1
        else:
            Use_Smin_x = min(sig_x1 , sig_x2)

        if sig_y1 * sig_y2 >= 0:
            Use_Smin_y = sig_y2 if abs(sig_y1) > abs(sig_y2) else sig_y1
        else:
            Use_Smin_y = min(sig_y1 , sig_y2)

        shear_ratio_long = 1 if Use_Smax_x == 0 else Use_Smin_x / Use_Smax_x
        shear_ratio_trans = 1 if Use_Smax_y == 0 else Use_Smin_y/Use_Smax_y

        Max_vonMises_x = sig_x1 if abs(sig_x1) > abs(sig_x2) else sig_x2

        unstf_pl_data['sxsd'] = sxsd
        unstf_pl_data['sy1sd'] = sy1sd

        l1 = min(l/4, s/2)
        if l == 0:
            sig_trans_l1 = Use_Smax_y
        else:
            sig_trans_l1 = Use_Smax_y*(shear_ratio_trans+(1-shear_ratio_trans)*(l-l1)/l)

        trans_stress_used = sysd = 0.75*Use_Smax_y if abs(0.75*Use_Smax_y) > abs(sig_trans_l1) else sig_trans_l1
        unstf_pl_data['sysd'] = sysd
        #Pnt. 5  Lateral loaded plates

        sjsd =math.sqrt(math.pow(Max_vonMises_x,2) + math.pow(sysd,2)-Max_vonMises_x*sysd+3*math.pow(tsd,2))

        uf_sjsd = sjsd/fy*gammaM
        unstf_pl_data['UF sjsd'] = uf_sjsd

        #psi_x =max([0,(1-math.pow(sjsd/fy,2))/math.sqrt(1-3/4*math.pow(sysd/fy,2)-3*math.pow(tsd/fy,2))])
        psi_x =max([0,(1-math.pow(sjsd/fy,2))/math.sqrt(1-3/4*math.pow(sysd/fy,2)-3*math.pow(tsd/fy,2))]) \
            if 1-3/4*math.pow(sysd/fy,2)-3*math.pow(tsd/fy,2) > 0 else 0
        psi_x_chk = (1-3/4*math.pow(sy1sd/fy,2)-3*math.pow(tsd/fy,2))>0

        psi_y = max([0,(1-math.pow(sjsd/fy,2))/math.sqrt(1-3/4*math.pow(sxsd/fy,2)-3*math.pow(tsd/fy,2))]) \
            if 1-3/4*math.pow(sxsd/fy,2)-3*math.pow(tsd/fy,2) > 0 else 0
        psi_y_chk = (1 - 3 / 4 * math.pow(sxsd / fy, 2) - 3 * math.pow(tsd / fy, 2)) > 0

        if gammaM * s * l == 0:
            Psd_max_press = 0
        else:
            if all([psi_x_chk, psi_y_chk]):
                Psd_max_press = (4 * fy / gammaM * math.pow(t / s,2) * (psi_y + math.pow(s / l, 2) * psi_x))
            else:
                Psd_max_press = -1

        if Psd_max_press == 0:
            uf_lat_load_pl_press = 0
        else:
            uf_lat_load_pl_press = 9 if Psd_max_press < 0 else abs(psd/Psd_max_press)

        unstf_pl_data['UF Pnt. 5  Lateral loaded plates'] = uf_lat_load_pl_press
        #6.2 & 6.6 Longitudinal stress
        if shear_ratio_long <= -2:
            ksig = "Unknown"
        elif 0 <= shear_ratio_long <= 1:
            ksig = 8.2 / (1.05 + shear_ratio_long)
        elif -1 <= shear_ratio_long < 0:
            ksig = 7.81 - 6.29 * shear_ratio_long + 9.78 * math.pow(shear_ratio_long, 2)
        elif -2 < shear_ratio_long < -1:
            ksig = 5.98 * math.pow(1 - shear_ratio_long, 2)

        #print(sxsd, sy1sd, tsd, sjsd, uf_lat_load_pl, psi_x, psi_y, uf_lat_load_pl_press, psd, Psd_max_press,ksig)

        if t*E == 0:
            alpha_p = 0
        elif ksig == "Unknown":
            alpha_p = 1.05*s/t*math.sqrt(fy/E)
        else:
            alpha_p = s/t/(28.4*math.sqrt(ksig*235/fy))

        if alpha_p <= 0:
            Cx = 0
        elif alpha_p <= 0.673:
            Cx = 1
        else:
            Cx = (alpha_p-0.055*(3+max([-2,shear_ratio_long])))/math.pow(alpha_p, 2)
            Cx = min(1, max(0, Cx))

        sxRd = Cx*fy/gammaM if not all([sig_x1<0, sig_x2<0]) else 1*fy/gammaM # Corrected 07.08.2023, issue 126

        uf_unstf_pl_long_stress = 0 if sxRd == 0 else abs(sxsd/sxRd)
        unstf_pl_data['UF Longitudinal stress'] = uf_unstf_pl_long_stress
        #print(uf_unstf_pl_long_stress)

        #6.3 & 6.8 Transverse stresses:
        ha = 0 if t == 0 else max([0,0.05*s/t-0.75])
        kp_1_for_Psd = 0 if s == 0 else 2*math.pow(t/s,2)*fy
        kp_used = 1-ha*(psd/fy-2*math.pow(t/s,2)) if psd>kp_1_for_Psd else 1

        alpha_c = 0 if t*E == 0 else 1.1*s/t*math.sqrt(fy/E)
        mu = 0.21*(alpha_c-0.2)

        if alpha_c <= 0.2:
            kappa = 1
        elif 0.2 < alpha_c < 2:
            kappa = 0 if alpha_c == 0 else 1/(2*math.pow(alpha_c,2))*(1+mu+math.pow(alpha_c,2)-
                                                                      math.sqrt(math.pow(1+mu+math.pow(alpha_c,2),2)-
                                                                                4*math.pow(alpha_c,2)))
        elif alpha_c >= 2:
            kappa = 0 if alpha_c == 0 else 1/(2*math.pow(alpha_c,2))+0.07


        syR = 0 if l*fy == 0 else (1.3*t/l*math.sqrt(E/fy)+kappa*(1-1.3*t/l*math.sqrt(E/fy)))*fy*kp_used
        syRd = syR if not all([sig_y1<0, sig_y2<0]) else fy
        syRd = syRd/gammaM
        uf_unstf_pl_trans_stress = 0 if syRd == 0 else abs(sysd)/syRd
        #print(uf_unstf_pl_trans_stress)
        unstf_pl_data['syR'] = syR
        unstf_pl_data['syRd'] = syRd
        unstf_pl_data['UF transverse stresses'] = uf_unstf_pl_trans_stress
        #6.4  Shear stress
        if l >= s:
            kl = 0 if l == 0 else 5.34+4*math.pow(s/l,2)
        else:
            kl = 0 if l == 0 else 5.34*math.pow(s/l,2)+4
        unstf_pl_data['kl'] = kl
        alpha_w = 0 if t*E*kl == 0 else 0.795*s/t*math.sqrt(fy/E/kl)
        if alpha_w <= 0.8:
            Ctau = 1
        elif 0.8 < alpha_w < 1.25:
            Ctau = 1-0.675*(alpha_w-0.8)
        else:
            Ctau = 0 if alpha_w == 0 else 0.9/alpha_w

        tauRd = Ctau*fy/gammaM/math.sqrt(3)
        uf_unstf_pl_shear_stress = 0 if tauRd == 0 else tsd/tauRd
        unstf_pl_data['UF Shear stresses'] = uf_unstf_pl_shear_stress
        #print(uf_unstf_pl_shear_stress)

        #6.5  Combined stresses

        if alpha_w <= 0.8:
            Ctaue = 1
        elif 0.8 < alpha_w < 1.25:
            Ctaue = 1-0.8*(alpha_w-0.8)
        else:
            Ctaue = 0 if alpha_w == 0 else 1/math.pow(alpha_w,2)

        tauRd_comb = Ctaue*fy/gammaM/math.sqrt(3)
        tauRd_comb = tauRd if sysd>0 else tauRd

        if s/t <= 120:
            ci = 0 if t == 0 else 1-s/120/t
        elif s/t > 120:
            ci  = 0
        else:
            ci = 1

        sxRd_comb = fy/gammaM if all([sig_x1<0, sig_x2<0]) else sxRd
        syRd_comb = syRd

        sxsd_div_sxrd = 0 if sxRd_comb == 0 else sxsd/sxRd_comb
        sysd_div_syrd = 0 if syRd_comb == 0 else sysd / syRd_comb
        tausd_div_taurd = 0 if tauRd_comb == 0 else tsd/tauRd_comb

        comb_req = math.pow(sxsd_div_sxrd, 2)+math.pow(sysd_div_syrd, 2)-ci*sxsd_div_sxrd*sysd_div_syrd+\
                   math.pow(tausd_div_taurd, 2)
        uf_unstf_pl_comb_stress = comb_req
        unstf_pl_data['UF Combined stresses'] = uf_unstf_pl_comb_stress

        return unstf_pl_data

    def stiffened_panel(self, unstf_pl_data = None, optimizing = False):
        E = self._E / 1e6
        v = self._v
        G = E/(2*(1+v))
        fy = self._mat_yield/1e6
        gammaM = self._Plate.mat_factor
        t = self._Plate.t
        s = self._Plate.spacing
        l = self._Plate.span * 1000

        sig_x1 = self._Plate.sigma_x1 * self._stress_load_factor
        sig_x2 = self._Plate.sigma_x2 * self._stress_load_factor

        sig_y1 = self._Plate.sigma_y1 * self._stress_load_factor
        sig_y2 = self._Plate.sigma_y2 * self._stress_load_factor

        Lg = self._Plate.girder_lg*1000

        stf_pl_data = dict()

        sxsd = unstf_pl_data['sxsd']
        #sxsd = 0 if self._method == 2 else unstf_pl_data['sxsd']
        sy1sd = 0 if self.get_method() == 2 else unstf_pl_data['sy1sd']

        sysd = 0 if self.get_method() == 2 else unstf_pl_data['sysd']
        tsd = abs(self._Plate.tau_xy * self._stress_load_factor)
        psd = self._lat_press * self._lat_load_factor
        psd_min_adj = psd if self._min_lat_press_adj_span is None else\
            self._min_lat_press_adj_span*self._lat_load_factor
        shear_ratio_long = 1
        shear_ratio_trans = 1

        #Pnt.7:  Buckling of stiffened plates
        Vsd = psd*s*l/2
        tw_req = Vsd*gammaM*math.sqrt(3)/(fy*self._Stiffener.hw)
        Anet = (self._Stiffener.hw + self._Stiffener.tf) * self._Stiffener.tw# + self._Stiffener.b*self._Stiffener.tf
        Vrd = Anet*fy/(gammaM*math.sqrt(3))
        Vsd_div_Vrd = Vsd/Vrd

        stf_pl_data['UF Shear force'] = Vsd_div_Vrd
        if optimizing and Vsd_div_Vrd > 1:
            return ['UF Shear force', Vsd_div_Vrd]

        # 7.2  Forces in idealised stiffened plate
        Iy = Is = self._Stiffener.get_moment_of_intertia()*1000**4

        stf_pl_data['Is'] = Is
        kc = 0 if t*s == 0 else 2*(1+math.sqrt(1+10.9*Is/(math.pow(t,3)*s)))
        mc = 13.3 if self._stf_end_support == 'Continuous' else 8.9

        # 7.3 Effective plate width
        syR = unstf_pl_data['syR']

        cys_tension_term = 4-3*math.pow(sysd/fy, 2)
        Cys = 0.5*((0 if cys_tension_term <= 0 else math.sqrt(cys_tension_term))+sysd/fy)

        alphap = 0 if t*E == 0 else 0.525 * (s / t) * math.sqrt(fy / E)  # reduced plate slenderness, checked not calculated with ex
        Cxs = (alphap - 0.22) / math.pow(alphap, 2) if alphap > 0.673 else 1
        stf_pl_data['alphap'] = alphap
        stf_pl_data['Cxs'] = Cxs
        if sysd < 0:
            Cys = max(0, min(Cys, 1))
        else:
            if s / t <= 120:
                ci = 0 if t == 0 else 1-s / 120 / t
            else:
                ci = 0

            cys_chk = 1 - math.pow(sysd / syR, 2) + ci * ((sxsd * sysd) / (Cxs * fy * syR))
            Cys =0 if cys_chk < 0 else math.sqrt(cys_chk)

        stf_pl_data['Cys_comp'] = Cys

        se_div_s = Cxs * Cys
        se = s * se_div_s

        zp = self._Stiffener.get_cross_section_centroid_with_effective_plate(se=se/1000)*1000 - t / 2  # ch7.5.1 page 19
        zt = (self._Stiffener.hw+self._Stiffener.tf) - zp + t/2

        Iy = self._Stiffener.get_moment_of_intertia(efficent_se=se/1000)*1000**4

        Weff = 0.0001 if zt == 0 else Iy/zt
        Co= 0 if kc*E*t*s == 0 else Weff*fy*mc/(kc*E*math.pow(t,2)*s)
        Po = 0 if all([sig_y1 < 0, sig_y2 < 0]) else (0.6+0.4*shear_ratio_trans)*Co*sy1sd \
            if shear_ratio_trans >-1.5 else 0

        qsd_press = (psd+abs(Po))*s
        qsd_opposite = abs(Po)*s if psd < Po else 0

        '''
        1	Overpressure on Stiffener Side
        2	Overpressure on Plate Side
        3	Overpr. may occur on both sides
        '''

        qsd_plate_side = qsd_opposite if self._overpressure_side == 'stiffener side' else qsd_press
        qsd_stf_side = qsd_opposite if self._overpressure_side == 'plate side' else qsd_press
        kl = unstf_pl_data['kl']

        tcrl = 0 if s == 0 else kl*0.904*E*math.pow(t/s,2)

        if l<= Lg:
            kg = 0 if Lg == 0 else 5.34+4*math.pow(l/Lg,2)
        else:
            kg = 0 if Lg == 0 else 5.34*math.pow(l / Lg, 2)+4

        tcrg = 0 if l == 0 else kg*0.904*E*math.pow(t/l,2)

        if self._tension_field_action == 'allowed' and tsd>(tcrl/gammaM):
            ttf = tsd-tcrg
        else:
            ttf = 0

        As = self._Stiffener.tw*self._Stiffener.hw + self._Stiffener.b*self._Stiffener.tf

        NSd = sxsd*(As+s*t)+ttf*s*t

        #7.4  Resistance of plate between stiffeners
        ksp = math.sqrt(1-3*math.pow(tsd/fy,2)) if tsd < (fy/math.sqrt(3)) else 0
        syrd_unstf = unstf_pl_data['syRd'] * ksp
        tsd_7_4 = fy/(math.sqrt(3)*gammaM)
        uf_stf_panel_res_bet_plate = max([sysd/syrd_unstf if all([syrd_unstf >0, sysd > 0]) else 0, tsd/tsd_7_4])
        stf_pl_data['UF Plate resistance'] = uf_stf_panel_res_bet_plate
        if optimizing and uf_stf_panel_res_bet_plate > 1:
            return ['UF Plate resistance', uf_stf_panel_res_bet_plate]
        #7.5  Characteristic buckling strength of stiffeners

        fEpx = 0 if s == 0 else 3.62*E*math.pow(t/s,2) # eq 7.42, checked, ok
        fEpy = 0 if s == 0 else 0.9*E*math.pow(t/s,2) # eq 7.43, checked, ok
        fEpt = 0 if s == 0 else 5.0*E*math.pow(t/s,2) # eq 7.44, checked, ok
        c = 0 if l == 0 else 2-(s/l) # eq 7.41, checked, ok

        sjSd = math.sqrt(
            math.pow(max([sxsd, 0]), 2) + math.pow(max([sysd, 0]), 2) - max([sxsd, 0]) * max([sysd, 0]) +
            3 * math.pow(tsd, 2))  # eq 7.38, ok

        if sjSd == 0:
            alphae = 0
        else:
            alphae = math.sqrt((fy/sjSd) * math.pow(math.pow(max([sxsd, 0])/fEpx, c)+
                                                    math.pow(max([sysd, 0])/fEpy, c)+
                                                    math.pow(abs(tsd)/fEpt, c), 1/c)) # eq 7.40



        fep = fy / math.sqrt(1+math.pow(alphae,4)) # eq 7.39
        eta = min(sjSd/fep, 1) # eq. 7.377

        C = 0 if self._Stiffener.tw == 0 else (self._Stiffener.hw / s) * math.pow(t / self._Stiffener.tw, 3) * \
                                              math.sqrt((1 - eta)) # e 7.36, checked ok

        beta = (3*C+0.2)/(C+0.2) # eq 7.35
        It = self._Stiffener.get_torsional_moment_venant(efficient_flange=False)
        Ipo = self._Stiffener.get_polar_moment()
        Iz = self._Stiffener.get_Iz_moment_of_inertia()

        def red_prop():
            tw_red =max(0,self._Stiffener.tw*(1-Vsd_div_Vrd))
            Atot_red  = As+se*t-self._Stiffener.hw*(self._Stiffener.tw - tw_red )
            It_red  = self._Stiffener.get_torsional_moment_venant(reduced_tw=tw_red, efficient_flange=False)
            Ipo_red  = self._Stiffener.get_polar_moment(reduced_tw=tw_red )
            #Iz = self._Stiffener.get_Iz_moment_of_inertia(reduced_tw=tw)
            #Iz_red = self._Stiffener.get_moment_of_intertia(efficent_se=se/1000, reduced_tw=tw_red)
            Iy_red = self._Stiffener.get_moment_of_intertia(efficent_se=se / 1000, reduced_tw=tw_red) * 1000 ** 4
            zp_red  = self._Stiffener.get_cross_section_centroid_with_effective_plate(se / 1000, reduced_tw=tw_red ) \
                      * 1000 - t / 2  # ch7.5.1 page 19
            zt_red  = (self._Stiffener.hw + self._Stiffener.tf) - zp_red + t/2  # ch 7.5.1 page 19
            Wes_red  = 0.0001 if zt_red == 0 else Iy_red/zt_red
            Wep_red  = 0.0001 if zp_red == 0 else Iy_red/zp_red
            return {'tw':tw_red , 'Atot': Atot_red , 'It': It_red , 'Ipo': Ipo_red , 'zp': zp_red ,
                    'zt': zt_red , 'Wes': Wes_red , 'Wep': Wep_red, 'Iy': Iy_red}

        hs = self._Stiffener.hw / 2 if self._Stiffener.get_stiffener_type() == 'FB' else \
            self._Stiffener.hw + self._Stiffener.tf / 2

        def lt_params(lT):

            if Ipo*lT>0:
                fET = beta*G*It/Ipo+math.pow(math.pi,2)*E*math.pow(hs,2)*Iz/(Ipo*math.pow(lT,2)) #NOTE, beta was missed from above, added. 23.08.2022
            else:
                fET = 0.001
            alphaT = 0 if fET == 0 else math.sqrt(fy/fET)
            mu = 0.35*(alphaT-0.6)
            fT_div_fy = (1+mu+math.pow(alphaT,2)-math.sqrt(math.pow(1+mu+math.pow(alphaT,2),2)-
                                                           4*math.pow(alphaT,2)))/(2*math.pow(alphaT,2))
            fT = fy*fT_div_fy if alphaT > 0.6 else fy
            #print(fET, alphaT, mu, fT)
            return {'fEt': fET, 'alphaT': alphaT, 'mu': mu, 'fT_div_fy': fT_div_fy, 'fT': fT}


        zp = self._Stiffener.get_cross_section_centroid_with_effective_plate(se/1000)*1000 - t / 2  # ch7.5.1 page 19
        zt  = (self._Stiffener.hw + self._Stiffener.tf) - zp + t/2  # ch 7.5.1 page 19
        fr = fy

        if Vsd_div_Vrd < 0.5:
            Wes = 0.0001 if zt == 0 else Iy/zt
            Wep = 0.0001 if zp == 0 else Iy/zp
            Ae = As + se * t
        else:
            red_param = red_prop()
            Wes = red_param['Wes']
            Wep = red_param['Wep']
            Ae = red_param['Atot']

        Wmin = min([Wes, Wep])
        pf = 0.0001 if l*s*gammaM == 0 else 12*Wmin*fy/(math.pow(l,2)*s*gammaM)

        if self._buckling_length_factor_stf is None:
            if self._stf_end_support == 'Continuous':
                lk = l*(1-0.5*abs(psd_min_adj/pf))

            else:
                lk = l
        else:
            lk = self._buckling_length_factor_stf * l
        ie = 0.0001 if As+se*t == 0 else math.sqrt(Iy/(As+se*t))
        fE = 0.0001 if lk == 0 else math.pow(math.pi,2)*E*math.pow(ie/lk,2)


        fk_dict = dict()
        fr_dict = dict()
        #Plate side
        zp = zp
        fr = fy
        fr_dict['plate'] = fr
        alpha = math.sqrt(fr/fE)
        mu = 0 if ie == 0 else (0.34+0.08*zp/ie)*(alpha-0.2)
        fk_div_fr = (1+mu+math.pow(alpha,2)-math.sqrt(math.pow(1+mu+math.pow(alpha,2),2)-4*math.pow(alpha,2)))/(2*math.pow(alpha,2))
        fk = fk_div_fr*fr if alpha > 0.2 else fr
        fk_dict['plate'] = fk
        #Stiffener side

        for lT in [int(l if self._stf_dist_between_lateral_supp is None else self._stf_dist_between_lateral_supp),
                   int(0.4*l if self._stf_dist_between_lateral_supp is None else self._stf_dist_between_lateral_supp),
                   int(0.8*l if self._stf_dist_between_lateral_supp is None else self._stf_dist_between_lateral_supp)]:
            params = lt_params(lT)
            fr = params['fT'] if params['alphaT']>0.6 else fy
            fr_dict[lT] = fr
            alpha = math.sqrt(fr / fE)
            mu = 0 if ie == 0 else (0.34 + 0.08 * zt  / ie) * (alpha - 0.2)
            fk_div_fr = (1 + mu + math.pow(alpha, 2) - math.sqrt(
                math.pow(1 + mu + math.pow(alpha, 2), 2) - 4 * math.pow(alpha, 2))) / (2 * math.pow(alpha, 2))
            fk = fk_div_fr * fr if alpha > 0.2 else fr
            fk_dict[lT] = fk

        #7.7.3  Resistance parameters for stiffeners

        NRd = 0.0001 if gammaM == 0 else Ae * (fy / gammaM)  # eq7.65, checked ok
        NksRd = Ae * (fk_dict[int(l if self._stf_dist_between_lateral_supp is None else self._stf_dist_between_lateral_supp)] / gammaM) #eq7.66
        NkpRd = Ae * (fk_dict['plate'] / gammaM)  # checked ok

        Ms1Rd = Wes * (fr_dict[int(0.4*l if self._stf_dist_between_lateral_supp is None else
                                   self._stf_dist_between_lateral_supp)] / gammaM)  # ok
        Ms2Rd = Wes * (fr_dict[int(0.8*l if self._stf_dist_between_lateral_supp is None else
                                   self._stf_dist_between_lateral_supp)] / gammaM)  # eq7.69 checked ok

        MstRd = Wes*(fy/gammaM) #eq7.70 checked ok
        MpRd = Wep*(fy/gammaM) #eq7.71 checked ok

        Ne = ((math.pow(math.pi,2))*E*Ae)/(math.pow(lk/ie,2))# eq7.72 , checked ok

        #7.6  Resistance of stiffened panels to shear stresses
        Ip = math.pow(t,3)*s/10.9
        tcrs = (36 * E / (s * t * math.pow(l, 2))) * ((Ip * math.pow(Is, 3)) ** 0.25)
        tRdy = fy/math.sqrt(3)/gammaM
        tRdl = tcrl/gammaM
        tRds = tcrs/gammaM
        tRd = min([tRdy,tRdl,tRds])


        u = 0 if all([tsd>(tcrl/gammaM), self._tension_field_action == 'allowed']) else math.pow(tsd/tRd, 2)
        zstar = zp
        if self._stf_end_support != 'Continuous':
            #Lateral pressure on plate side:
            #7.7.2 Simple supported stiffener (sniped stiffeners)

            #Lateral pressure on plate side:
            stf_pl_data['UF Stiffener side'] = 0
            stf_pl_data['UF Plate side'] = 0
            uf_7_58 = NSd/NksRd-2*NSd/NRd +((qsd_plate_side*math.pow(l,2)/8)+NSd*zstar)/(MstRd*(1-NSd/Ne))+u
            uf_7_59 = NSd/NkpRd+((qsd_plate_side*math.pow(l,2)/8)+NSd*zstar)/(MpRd*(1-NSd/Ne))+u
            uf_max_simp_pl = max([uf_7_58, uf_7_59])
            stf_pl_data['UF simply supported plate side'] = uf_max_simp_pl

            #Lateral pressure on stiffener side:

            uf_7_60 = NSd/NksRd+((qsd_stf_side*math.pow(l,2)/8)-NSd*zstar)/(Ms2Rd*(1-NSd/Ne))+u
            uf_7_61 = NSd/NkpRd-2*NSd/NRd+((qsd_stf_side*math.pow(l,2)/8)-NSd*zstar)/(MpRd*(1-NSd/Ne))+u

            test_qsd_l = qsd_stf_side*math.pow(l,2)/8 >= NSd*zstar
            uf_7_62 = NSd/NksRd-2*NSd/NRd+(NSd*zstar-(qsd_stf_side*math.pow(l,2)/8))/(MstRd*(1-NSd/Ne))+u
            uf_7_63 = NSd/NkpRd+(NSd*zstar-(qsd_stf_side*math.pow(l,2)/8))/(MpRd*(1-NSd/Ne))+u

            uf_max_simp_stf = max([0,uf_7_62,uf_7_63]) if not test_qsd_l else max([0,uf_7_60,uf_7_61])
            stf_pl_data['UF simply supported stf side'] = uf_max_simp_stf
        else:
            stf_pl_data['UF simply supported stf side'] = 0
            stf_pl_data['UF simply supported plate side'] = 0
            #7.7.1 Continuous stiffeners

            M1Sd_pl = abs(qsd_plate_side)*math.pow(l,2)/self._km3
            M2Sd_pl = abs(qsd_plate_side)*math.pow(l,2)/self._km2
            M1Sd_stf = abs(qsd_stf_side) * math.pow(l, 2) / self._km3
            M2Sd_stf = abs(qsd_stf_side) * math.pow(l, 2) / self._km2
            M1Sd_max = max([M1Sd_pl, M1Sd_stf])
            M2Sd_max = max([M2Sd_pl, M2Sd_stf])
            # Lateral pressure on plate side:
            #print(M1Sd_pl, M2Sd_pl, M1Sd_stf,M2Sd_stf, qsd_stf_side, qsd_plate_side)
            from scipy.optimize import minimize
            def iteration_min_uf_pl_side(x):
                eq7_50 = NSd/NksRd+(M1Sd_pl-NSd*x)/(Ms1Rd*(1-NSd/Ne))+u
                eq7_51 = NSd/NkpRd-2*NSd/NRd +(M1Sd_pl-NSd*x)/(MpRd*(1-NSd/Ne))+u
                eq7_52 = NSd/NksRd-2*NSd/NRd+(M2Sd_pl+NSd*x)/(MstRd*(1-NSd/Ne))+u
                eq7_53 = NSd/NkpRd+(M2Sd_pl+NSd*x)/(MpRd*(1-NSd/Ne))+u
                #print(zstar, eq7_50, eq7_51,eq7_52,eq7_53,max([eq7_50, eq7_51,eq7_52,eq7_53]))
                return max(eq7_50, eq7_51, eq7_52, eq7_53)
            res_iter_pl = minimize(iteration_min_uf_pl_side, 0, bounds=[[-zt+self._Stiffener.tf/2,zp]])

            if type(res_iter_pl.fun) == list:
                stf_pl_data['UF Plate side'] = res_iter_pl.fun[0]
            else:
                stf_pl_data['UF Plate side'] = res_iter_pl.fun

            # Lateral pressure   on stiffener side:

            # max_lfs = []
            # ufs = []

            def iteration_min_uf_stf_side(x):
                eq7_54 = NSd/NksRd-2*NSd/NRd +(M1Sd_stf+NSd*x)/(MstRd*(1-NSd/Ne))+u
                eq7_55 = NSd/NkpRd+(M1Sd_stf+NSd*x)/(MpRd*(1-NSd/Ne))+u
                eq7_56 = NSd/NksRd+(M2Sd_stf-NSd*x)/(Ms2Rd*(1-NSd/Ne))+u
                eq7_57 = NSd/NkpRd-2*NSd/NRd+(M2Sd_stf-NSd*x)/(MpRd*(1-NSd/Ne))+u
                return max(eq7_54, eq7_55, eq7_56, eq7_57)

            res_iter_stf = minimize(iteration_min_uf_stf_side, 0, bounds=[[-zt+self._Stiffener.tf/2,zp]])

            if type(res_iter_stf.fun) == list:
                stf_pl_data['UF Stiffener side'] = res_iter_stf.fun[0]
            else:
                stf_pl_data['UF Stiffener side'] = res_iter_stf.fun

        return stf_pl_data

    def girder(self, unstf_pl_data = None, stf_pl_data = None, optmizing = False):
        '''
        Buckling of girder.
        '''

        girder_data = dict()

        E = self._E / 1e6
        v = self._v
        G = E/(2*(1+v))
        fy = self._mat_yield/1e6
        gammaM = self._Plate.mat_factor
        t = self._Plate.t
        s = self._Plate.spacing
        l = self._Plate.span * 1000
        hw = self._Girder.hw

        tsd = abs(self._Plate.tau_xy * self._stress_load_factor)
        psd = self._lat_press


        sig_x1 = self._Plate.sigma_x1 * self._stress_load_factor
        sig_x2 = self._Plate.sigma_x2 * self._stress_load_factor

        sig_y1 = self._Plate.sigma_y1 * self._stress_load_factor
        sig_y2 = self._Plate.sigma_y2 * self._stress_load_factor

        sxsd = unstf_pl_data['sxsd']
        #sxsd = 0 if self._method == 2 else unstf_pl_data['sxsd']

        sy1sd = unstf_pl_data['sy1sd']

        #sysd = 0 if self.get_method() == 2 else unstf_pl_data['sysd']
        sysd = unstf_pl_data['sysd']
        tsd = abs(self._Plate.tau_xy * self._stress_load_factor)
        psd = self._lat_press * self._lat_load_factor
        psd_min_adj = psd if self._min_lat_press_adj_span is None else\
            self._min_lat_press_adj_span*self._lat_load_factor

        Lg = self._Plate.girder_lg*1000

        Ltg = Lg if self._girder_dist_between_lateral_supp == None else self._girder_dist_between_lateral_supp
        Lp = 0 if self._panel_length_Lp is None else self._panel_length_Lp

        #Pnt.8:  Buckling of Girders
        #7.8  Check for shear force
        Vsd = psd*l*Lg/2

        tw_req = Vsd*gammaM*math.sqrt(3)/(fy*self._Girder.hw)
        Anet = self._Girder.hw * self._Girder.tw + self._Girder.tw*self._Girder.tf
        Vrd = Anet*fy/(gammaM*math.sqrt(3))

        Vsd_div_Vrd = Vsd/Vrd
        girder_data['UF shear force'] = Vsd_div_Vrd
        if optmizing and Vsd_div_Vrd > 1:
            return ['UF shear force', Vsd_div_Vrd]

        CHK_account_for_interaction = Vsd < 0.5*Vrd

        #8.2  Girder forces
        As = self._Stiffener.tw*self._Stiffener.hw + self._Stiffener.b*self._Stiffener.tf
        Ag = self._Girder.tw*self._Girder.hw + self._Girder.b*self._Girder.tf


        #sysd = 0 if self.get_method() == 2 else unstf_pl_data['sysd']
        NySd = sysd*(Ag+l*t)

        Is = stf_pl_data['Is']

        tcel = 18*E/(t*math.pow(l,2))*math.pow(t*Is/s, 0.75)
        tceg = 0 if Lp == 0 else tcel*math.pow(l,2)/math.pow(Lp,2)

        alpha_t1 = 0 if Lp == 0 else math.sqrt(0.6*fy/tceg)
        alpha_t2 = math.sqrt(0.6*fy/tcel)

        tcrg = 0.6*fy/math.pow(alpha_t1,2) if alpha_t1 > 1 else 0.6*fy
        tcrl = 0.6*fy/math.pow(alpha_t2,2) if alpha_t2 > 1 else 0.6*fy

        tcrg = tcrg if self._stf_end_support == 'Continuous' else 0

        #8.4 Effective width of girders

        #Method 1:
        alphap = stf_pl_data['alphap']
        Cxs = stf_pl_data['Cxs']
        fkx = Cxs*fy
        sxsd_compression = max(sxsd, 0)
        CxG = math.sqrt(max(0, 1-math.pow(sxsd_compression/fkx, 2))) if fkx > 0 else 0
        if 4-math.pow(Lg/l,2) != 0:
            CyG_tens = 1 if Lg > 2*l else Lg/(l*math.sqrt(4-math.pow(Lg/l,2)))
        else:
            CyG_tens = 1
        CyG_comp  = 0 if l*alphap == 0 else stf_pl_data['Cys_comp']
        CyG = min([1,CyG_tens]) if sy1sd<0 else min([1, CyG_comp])
        CtG = math.sqrt(1-3*math.pow(tsd/fy,2)) if tsd<fy/math.sqrt(3) else 0
        le_div_l_method1 = CxG*CyG*CtG
        le_method1 = l*le_div_l_method1

        lim_sniped_or_cont = 0.3*Lg if self._girder_end_support == 'Continuous' else 0.4*Lg
        tot_min_lim = min([le_method1, lim_sniped_or_cont])


        #Method 2:
        CxG = math.sqrt(max(0, 1-math.pow(sxsd_compression/fy, 2))) if fy > 0 else 0
        alphaG = 0 if E*t == 0 else 0.525*l/t*math.sqrt(fy/E)
        CyG = (alphaG-0.22)/math.pow(alphaG,2) if alphaG>0.673 else 1
        CtG = math.sqrt(1-3*math.pow(tsd/fy,2)) if tsd<fy/math.sqrt(3) else 0
        le_div_l_method2  = CxG*CyG*CtG
        le_method2 = le_div_l_method2*l

        eff_width_sec_mod = tot_min_lim if self.get_method() == 1 else le_method2
        eff_width_other_calc = le_method1 if self.get_method() == 1 else le_method2

        le = eff_width_other_calc

        AtotG = Ag + le * t

        Iy = self._Girder.get_moment_of_intertia(efficent_se=le / 1000) * 1000 ** 4
        zp = self._Girder.get_cross_section_centroid_with_effective_plate(le / 1000) * 1000 - t / 2  # ch7.5.1 page 19
        zt = (t / 2 + self._Girder.hw + self._Girder.tf) - zp  # ch 7.5.1 page 19
#
        def red_prop():
            twG =max(0,self._Girder.tw*(1-Vsd_div_Vrd))

            le = eff_width_other_calc

            AtotG = Ag+le*t-self._Girder.hw*(self._Girder.tw - twG)
            Ipo = self._Girder.get_polar_moment(reduced_tw=twG)
            IzG = self._Girder.get_Iz_moment_of_inertia(reduced_tw=twG)
            Iy = self._Girder.get_moment_of_intertia(efficent_se=le/1000, reduced_tw=twG)*1000**4

            zp = self._Girder.get_cross_section_centroid_with_effective_plate(le / 1000, reduced_tw=twG) * 1000 - t / 2  # ch7.5.1 page 19
            zt = (t / 2 + self._Girder.hw + self._Girder.tf) - zp  # ch 7.5.1 page 19
            WeG = 0.0001 if zt == 0 else Iy / zt
            Wep = 0.0001 if zp == 0 else Iy / zp
            #print('In reduced', 'zp',zp,'zt',zt,'WeG',WeG,'Wep',Wep, 'Iy', Iy)
            return {'tw':twG, 'Atot': AtotG, 'Ipo': Ipo, 'IzG': IzG, 'zp': zp, 'zt': zt, 'WeG': WeG, 'Wep': Wep}

        if Vsd_div_Vrd < 0.5:
            WeG = 0.0001 if zt == 0 else Iy/zt
            Wep = 0.0001 if zp == 0 else Iy/zp
            AeG = Ag+eff_width_other_calc*t
        else:
            red_param = red_prop()
            WeG = red_param['WeG']
            Wep = red_param['Wep']
            AeG = red_param['Atot']

        # #from: 7.7.3  Resistance parameters for stiffeners
        Wmin = min([WeG, Wep])
        pf = 0.0001 if l*s*gammaM == 0 else 12*Wmin*fy/(math.pow(l,2)*s*gammaM)

        lk = Lg
        LGk = lk if self._buckling_length_factor_girder is None else lk*self._buckling_length_factor_girder

        #ie = 0 if Vsd_div_Vrd<0.5 else math.sqrt(Iy/AtotG)
        ie = math.sqrt(Iy / AtotG)
        fE = 0 if LGk == 0 else math.pow(math.pi,2)*E*math.pow(ie/LGk,2)

        # 8.2  Girder forces, cont
        alphaG = 0 if fE == 0 else math.sqrt(fy/fE)
        Q = 0 if alphaG-0.2<0 else min([1, alphaG-0.2])
        C_for_tsd_trg = Q*(7-5*math.pow(s/l,2))*math.pow((tsd-tcrg)/tcrl,2)
        C = C_for_tsd_trg if tsd>tcrg else 0
        p0lim = 0.02*(t+As/s)/l*(sxsd+C*tsd)
        p0calc = 0 if s*self._Girder.hw*Lg*E*l==0 else 0.4*(t+As/s)/(self._Girder.hw*(1-s/Lg))*fy/E*math.pow(Lg/l,2)\
                                                       *(sxsd+C*tsd)
        p0_compression = max([p0lim,p0calc])
        p0_tension = 0 if s*Lg*self._Girder.hw*E*l==0 else 0.4*(t+As/s)/(self._Girder.hw*(l-s/Lg))*gammaM/E\
                                                          *math.pow(Lg/l,2)*(C*tsd)
        p0 = p0_tension if sxsd<0 else p0_compression

        qSd_pressure = (psd+p0_tension)*l if sxsd<0 else (psd+p0_compression)*l
        qsd_oppsite = p0*l if psd<p0 else 0
        qSd_plate_side = qsd_oppsite if self._overpressure_side == 'stiffener side' else qSd_pressure
        qSd_girder_side = qsd_oppsite if self._overpressure_side == 'plate side' else qSd_pressure

        #8.5  Torsional buckling of girders
        Af = self._Girder.tf*self._Girder.b
        Aw = self._Girder.hw*self._Girder.tw
        Iz = self._Girder.get_Iz_moment_of_inertia()

        b = max([self._Girder.b, self._Girder.tw])
        C = 0.55 if self._Girder.get_stiffener_type() in ['T', 'FB'] else 1.1
        LGT0 = b*C*math.sqrt(E*Af/(fy*(Af+Aw/3))) #TODO can add a automatic check/message if torsional buckling shall be considered
        girder_data['Torsional buckling'] = 'Torsional buckling to be considered' if Ltg >LGT0 else \
            "Torsional buckling need not to be considered"
        def lt_params(LTG):
            fETG = math.pow(math.pi, 2)*E*Iz/((Af+Aw/3)*math.pow(LTG, 2))
            alphaTG = math.sqrt(fy/fETG)
            mu = 0.35*(alphaTG-0.6)
            fT_div_fy = (1 + mu + math.pow(alphaTG, 2) - math.sqrt(
                math.pow(1 + mu + math.pow(alphaTG, 2), 2) - 4 * math.pow(alphaTG, 2))) / (2 * math.pow(alphaTG, 2))
            fT = fT_div_fy*fy if alphaTG>0.6 else fy
            return {'fETG': fETG, 'alphaT': alphaTG, 'mu': mu, 'fT_div_fy': fT_div_fy, 'fT': fT}

        fk_dict = dict()
        fr_dict = dict()
        for lT in ['plate', Ltg, 0.4*Lg, 0.8*Lg]:
            if lT != 'plate':
                params = lt_params(lT)
                fr = params['fT'] if params['alphaT']>0.6 else fy
                alpha = math.sqrt(fr / fE)
                mu = 0 if ie == 0 else (0.34 + 0.08 * zt / ie) * (alpha - 0.2)
            else:
                fr = fy
                alpha = math.sqrt(fr / fE)
                mu = 0 if ie == 0 else (0.34 + 0.08 * zp / ie) * (alpha - 0.2)
            fr_dict[lT] = fr
            fk_div_fr = (1 + mu + math.pow(alpha, 2) - math.sqrt(
                math.pow(1 + mu + math.pow(alpha, 2), 2) - 4 * math.pow(alpha, 2))) / (2 * math.pow(alpha, 2))
            fk = fk_div_fr * fr if alpha > 0.2 else fr
            fk_dict[lT] = fk

        # #7.7.3  Resistance parameters for stiffeners

        NRd = 0.0001 if gammaM == 0 else AeG * (fy / gammaM)  # eq7.65, checked ok

        NksRd = AeG * (fk_dict[Ltg] / gammaM) #eq7.66
        NkpRd = AeG * (fk_dict['plate'] / gammaM)  # checked ok
        MsRd = WeG*fr_dict[Ltg]/gammaM
        Ms1Rd = WeG * (fr_dict[0.4*Lg] / gammaM)  # ok
        Ms2Rd = WeG * (fr_dict[0.8*Lg] / gammaM)  # eq7.69 checked ok

        MstRd = WeG*(fy/gammaM) #eq7.70 checked ok
        MpRd = Wep*(fy/gammaM) #eq7.71 checked ok

        NE = ((math.pow(math.pi,2))*E*AeG)/(math.pow(LGk/ie,2))# eq7.72 , checked ok
        # print(fr_dict)
        # print(fk_dict)
        # print('WeG', WeG, 'Wep', Wep)
        # print('NRd',NRd, 'NksRd',NksRd, 'NkpRd',NkpRd,'MsRd', MsRd,'MstRd', MstRd, 'Ms1Rd', Ms1Rd, 'Ms2Rd', Ms2Rd, 'MstRd', MstRd, 'MpRd', MpRd, 'Ne', Ne)

        #7.7  Interaction formulas for axial compression and lateral pressure
        #7.7.2 Simple supported girder (sniped girders)
        if self._girder_end_support != 'Continuous':
            u = 0
            zstar = zp
            girder_data['UF Cont. plate side'] = 0
            girder_data['UF Cont. girder side'] = 0

            # Lateral pressure on plate side:
            uf_7_58 = NySd/NksRd-2*NySd/NRd +((qSd_plate_side*math.pow(Lg, 2)/8)+NySd*zstar)/(MstRd*(1-NySd/NE))+u
            uf_7_59 = NySd/NkpRd+((qSd_plate_side*math.pow(Lg, 2)/8)+NySd*zstar)/(MpRd*(1-NySd/NE))+u

            max_uf_simp_plate = max([0,uf_7_58,uf_7_59])
            girder_data['UF Simplified plate side'] = max_uf_simp_plate

            #Lateral pressure on girder side:
            uf_7_60 = NySd/NksRd+((qSd_girder_side*math.pow(Lg, 2)/8)-NySd*zstar)/(Ms2Rd*(1-NySd/NE))+u
            uf_7_61 = NySd/NkpRd-2*NySd/NRd+((qSd_girder_side*math.pow(Lg, 2)/8)-NySd*zstar)/(MpRd*(1-NySd/NE))+u

            CHK_qSd_NSd = qSd_girder_side*math.pow(Lg, 2)/8 < NySd*zstar

            uf_7_62 = NySd/NksRd-2*NySd/NRd+(NySd*zstar-(qSd_girder_side*math.pow(Lg, 2)/8))/(MstRd*(1-NySd/NE))+u
            uf_7_63 = NySd/NkpRd+(NySd*zstar-(qSd_girder_side*math.pow(Lg, 2)/8))/(MpRd*(1-NySd/NE))+u

            max_uf_simp_stiffener = max([0,uf_7_60,uf_7_61]) if CHK_qSd_NSd else max([0,uf_7_60,uf_7_61, uf_7_62,uf_7_63])
            girder_data['UF Simplified girder side'] = max_uf_simp_stiffener
        else:
            u = 0
            girder_data['UF Simplified girder side'] = 0
            girder_data['UF Simplified plate side'] = 0
            #7.7.1 Continuous stiffeners
            M1Sd_pl = abs(qSd_plate_side)*math.pow(Lg, 2)/12
            M2Sd_pl = abs(qSd_plate_side)*math.pow(Lg, 2)/24

            M1Sd_stf = abs(qSd_girder_side)*math.pow(Lg, 2)/12
            M2Sd_stf = abs(qSd_girder_side)*math.pow(Lg, 2)/24
            # #Lateral pressure on plate side:
            def iter_plate(zstar):
                uf_7_48 = NySd/NksRd+(M1Sd_pl-NySd*zstar)/(Ms1Rd*(1-NySd/NE))+u
                uf_7_49 = NySd/NkpRd-2*NySd/NRd +(M1Sd_pl-NySd*zstar)/(MpRd*(1-NySd/NE))+u
                uf_7_50 = NySd/NksRd-2*NySd/NRd+(M2Sd_pl+NySd*zstar)/(MstRd*(1-NySd/NE))+u
                uf_7_51 = NySd/NkpRd+(M2Sd_pl+NySd*zstar)/(MpRd*(1-NySd/NE))+u
                return max([uf_7_48, uf_7_49, uf_7_50, uf_7_51])

            res_iter_pl = minimize(iter_plate, 0, bounds=[[-zt + self._Girder.tf / 2, zp]])

            if type(res_iter_pl.fun) == list:
                girder_data['UF Cont. plate side'] = res_iter_pl.fun[0]
            else:
                girder_data['UF Cont. plate side'] = res_iter_pl.fun
            #     Lateral pressure on girder side:
            def iter_girder(zstar):
                uf_7_52 = NySd/NksRd-2*NySd/NRd +(M1Sd_stf +NySd*zstar)/(MstRd*(1-NySd/NE))+u
                uf_7_53 = NySd/NkpRd+(M1Sd_stf +NySd*zstar)/(MpRd*(1-NySd/NE))+u
                uf_7_54 = NySd/NksRd+(M2Sd_stf-NySd*zstar)/(Ms2Rd*(1-NySd/NE))+u
                uf_7_55 = NySd/NkpRd-2*NySd/NRd+(M2Sd_stf-NySd*zstar)/(MpRd*(1-NySd/NE))+u
                return max([uf_7_52, uf_7_53 ,uf_7_54 ,uf_7_55])

            res_iter_girder = minimize(iter_girder, 0, bounds=[[-zt + self._Girder.tf / 2, zp]])

            if type(res_iter_girder.fun) == list:
                girder_data['UF Cont. girder side'] = res_iter_girder.fun[0]
            else:
                girder_data['UF Cont. girder side'] = res_iter_girder.fun

        return girder_data

    def local_buckling(self, optimizing = False):
        '''
        Checks for girders and stiffeners
        '''
        fy = self._mat_yield / 1e6
        if self._Stiffener is not None:
            max_web_stf = 42*self._Stiffener.tw*math.sqrt(235/fy) if self._Stiffener.get_stiffener_type() != 'FB' else 0
            max_flange_stf = (14 if self._fab_method_stiffener == 'welded' else 15) *  self._Stiffener.tf *math.sqrt(235/fy)
        else:
            max_web_stf = 0
            max_flange_stf = 0

        if self._Girder is not None:
            max_web_girder = 42*self._Girder.tw*math.sqrt(235/fy) if self._Girder.get_stiffener_type() != 'FB' else 0
            max_flange_girder = (14 if self._fab_method_girder == 'welded' else 15) *  self._Girder.tf *math.sqrt(235/fy)
        else:
            max_web_girder = 0
            max_flange_girder = 0

        return {'Stiffener': [max_web_stf, max_flange_stf], 'Girder': [max_web_girder, max_flange_girder]}

    # def get_tuple(self):
    #     ''' Return a tuple of the plate stiffener'''
    #     return (self.Plate.get_s(), self.Plate.get_pl_thk(), self.Stiffener.get_web_thk(), self.Stiffener.get_web_thk(),
    #             self.Stiffener.get_fl_w(), self.Stiffener.get_fl_thk(), self.Plate.span, self.Plate.girder_lg,
    #             self.Stiffener.get_stiffener_type())

    def get_one_line_string_mixed(self):
        ''' Returning a one line string. '''
        return 'pl_'+str(round(self.Plate.get_s(), 1))+'x'+str(round(self.Plate.get_pl_thk(),1))+' stf_'+\
               self.Stiffener.get_stiffener_type()+\
               str(round(self.Stiffener.get_web_h(),1))+'x'+str(round(self.Stiffener.get_web_thk(),1))+'+'\
               +str(round(self.Stiffener.get_fl_w(),1))+'x'+\
               str(round(self.Stiffener.get_fl_thk(),1))

    def get_extended_string_mixed(self):
        ''' Some more information returned. '''
        return 'span: '+str(round(self.Plate.get_span(),4))+' structure type: '+ self.Stiffener.get_structure_type() + ' stf. type: ' + \
               self.Stiffener.get_stiffener_type() + ' pressure side: ' + self.overpressure_side

