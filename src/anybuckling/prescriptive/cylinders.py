'''
Prescriptive buckling of cylinders and curved plates according to
DNV-RP-C202 / DNV-OS-C101.

Transferred from ANYstructure (anystruct/calc_structure.py).
Contains the shell geometry model (Shell) and the cylinder buckling
calculation (CylinderAndCurvedPlate).
'''
import math
import numpy as np

from .plates import CalcScantlings, Structure

class Shell():
    '''
    Small class to contain shell properties.
    '''
    def __init__(self, main_dict: dict = None):
        super(Shell, self).__init__()
        '''
                            shell_dict = {'plate_thk': [self._new_shell_thk.get() / 1000, 'm'],
                                  'radius': [self._new_shell_radius.get() / 1000, 'm'],
                                  'distance between rings, l': [self._new_shell_dist_rings.get() / 1000, 'm'],
                                  'length of shell, L': [self._new_shell_length.get() / 1000, 'm'],
                                  'tot cyl length, Lc': [self._new_shell_tot_length.get() / 1000, 'm'],
                                  'eff. buckling lenght factor': [self._new_shell_k_factor.get() / 1000, 'm'],
                                  'mat_yield': [self._new_shell_yield.get() *1e6, 'Pa']}
        '''
        if main_dict is None:
            self._thk = None
            self._mat_yield = None
            self._radius = None
            self._dist_between_rings = None
            self._length_of_shell = None
            self._tot_cyl_length = None
            self._k_factor = None
        else:
            self._thk = main_dict['plate_thk'][0]
            self._mat_yield = main_dict['mat_yield'][0]
            self._radius = main_dict['radius'][0]
            self._dist_between_rings = main_dict['distance between rings, l'][0]
            self._length_of_shell = main_dict['length of shell, L'][0]
            self._tot_cyl_length = main_dict['tot cyl length, Lc'][0]
            self._k_factor = main_dict['eff. buckling lenght factor'][0]

        # For conical
        self._cone_r1 = None if main_dict is None else main_dict.get('cone r1', [None, 'm'])[0]
        self._cone_r2 = None if main_dict is None else main_dict.get('cone r2', [None, 'm'])[0]
        self._cone_length = None if main_dict is None else main_dict.get('cone length, l', [None, 'm'])[0]
        self._cone_alpha = None if main_dict is None else main_dict.get('cone alpha', [None, 'deg'])[0]
    @property
    def Lc(self):
        return self._tot_cyl_length
    @Lc.setter
    def Lc(self, val):
        self._tot_cyl_length = val
    @property
    def thk(self):
        return self._thk
    @thk.setter
    def thk(self, val):
        self._thk = val
    @property
    def radius(self):
        return self._radius
    @radius.setter
    def radius(self, val):
        self._radius = val
    @property
    def dist_between_rings(self):
        return self._dist_between_rings
    @dist_between_rings.setter
    def dist_between_rings(self, val):
        self._dist_between_rings = val
    @property
    def length_of_shell(self):
        return self._length_of_shell
    @length_of_shell.setter
    def length_of_shell(self, val):
        self._length_of_shell = val
    @property
    def tot_cyl_length(self):
        return self._tot_cyl_length
    @tot_cyl_length.setter
    def tot_cyl_length(self, val):
        self._tot_cyl_length = val
    @property
    def k_factor(self):
        return self._k_factor
    @k_factor.setter
    def k_factor(self, val):
        self._k_factor = val
    @property
    def cone_r1(self):
        return self._cone_r1
    @cone_r1.setter
    def cone_r1(self, val):
        self._cone_r1 = val
    @property
    def cone_r2(self):
        return self._cone_r2
    @cone_r2.setter
    def cone_r2(self, val):
        self._cone_r2 = val
    @property
    def cone_length(self):
        return self._cone_length
    @cone_length.setter
    def cone_length(self, val):
        self._cone_length = val
    @property
    def cone_alpha(self):
        return self._cone_alpha
    @cone_alpha.setter
    def cone_alpha(self, val):
        self._cone_alpha = val
    def set_conical_geometry(self, r1, r2, length):
        self._cone_r1 = r1
        self._cone_r2 = r2
        self._cone_length = length
        self._cone_alpha = math.degrees(math.atan(abs(r2 - r1) / length)) if length else 0
        self._radius = self.cone_equivalent_radius()
        self._dist_between_rings = self.cone_equivalent_length()
        self._length_of_shell = self._dist_between_rings
        self._tot_cyl_length = self._dist_between_rings
    def cone_equivalent_radius(self):
        if None in [self._cone_r1, self._cone_r2, self._cone_alpha]:
            return self._radius
        cos_alpha = math.cos(math.radians(self._cone_alpha))
        return (self._cone_r1 + self._cone_r2) / (2 * cos_alpha) if cos_alpha else 0
    def cone_equivalent_length(self):
        if None in [self._cone_length, self._cone_alpha]:
            return self._dist_between_rings
        cos_alpha = math.cos(math.radians(self._cone_alpha))
        return self._cone_length / cos_alpha if cos_alpha else 0
    def get_Zl(self):
        L = self.tot_cyl_length*1000
        Zl = math.pow(L,2)*math.sqrt(1-math.pow(0.3,2))/(self._radius*1000 * self._thk*1000) if self._thk*self._radius else 0
        return Zl

    def get_effective_width_shell_plate(self):
        return 1.56*math.sqrt(self._radius * self._thk)/(1+12*self.thk/self._radius)

    def get_main_properties(self):
        main_data = {'plate_thk': [self._thk, 'm'],
                                  'radius': [self._radius, 'm'],
                                  'distance between rings, l': [self._dist_between_rings, 'm'],
                                  'length of shell, L': [self._length_of_shell, 'm'],
                                  'tot cyl length, Lc': [self._tot_cyl_length, 'm'],
                                  'eff. buckling lenght factor': [self._k_factor, 'm'],
                                  'mat_yield': [self._mat_yield, 'Pa'],
                                  'cone r1': [self._cone_r1, 'm'],
                                  'cone r2': [self._cone_r2, 'm'],
                                  'cone length, l': [self._cone_length, 'm'],
                                  'cone alpha': [self._cone_alpha, 'deg']}
        return main_data

    def set_main_properties(self, main_dict):

        self._thk = main_dict['plate_thk'][0]
        self._mat_yield = main_dict['mat_yield'][0]
        self._radius = main_dict['radius'][0]
        self._dist_between_rings = main_dict['distance between rings, l'][0]
        self._length_of_shell = main_dict['length of shell, L'][0]
        self._tot_cyl_length = main_dict['tot cyl length, Lc'][0]
        self._k_factor = main_dict['eff. buckling lenght factor'][0]
        self._cone_r1 = main_dict.get('cone r1', [None, 'm'])[0]
        self._cone_r2 = main_dict.get('cone r2', [None, 'm'])[0]
        self._cone_length = main_dict.get('cone length, l', [None, 'm'])[0]
        self._cone_alpha = main_dict.get('cone alpha', [None, 'deg'])[0]

class CylinderAndCurvedPlate():
    '''
    Buckling of cylinders and curved plates.
    Geomeries
    	Selections for: Type of Structure Geometry:
    geomeries = {1:'Unstiffened shell (Force input)', 2:'Unstiffened panel (Stress input)',
                    3:'Longitudinal Stiffened shell  (Force input)',
                    4:'Longitudinal Stiffened panel (Stress input)',
                    5:'Ring Stiffened shell (Force input)',
                    6:'Ring Stiffened panel (Stress input)',
                    7:'Orthogonally Stiffened shell (Force input)',
                    8:'Orthogonally Stiffened panel (Stress input)'}

    '''

    geomeries = {11:'Flat plate, stiffened',10: 'Flat plate, unstiffened', 12: 'Flat plate, stiffened with girder',
                 1:'Unstiffened shell (Force input)', 2:'Unstiffened panel (Stress input)',
                 3:'Longitudinal Stiffened shell  (Force input)', 4:'Longitudinal Stiffened panel (Stress input)',
                 5:'Ring Stiffened shell (Force input)', 6:'Ring Stiffened panel (Stress input)',
                 7:'Orthogonally Stiffened shell (Force input)', 8:'Orthogonally Stiffened panel (Stress input)',
                 9:'Unstiffened conical shell (Force input)'}
    geomeries_map = dict()
    for key, value in geomeries.items():
        geomeries_map[value] = key

    geomeries_map_no_input_spec = dict()
    for key, value in geomeries.items():
        this_str = value.replace(' (Stress input)', '')
        this_str = this_str.replace(' (Force input)', '')
        geomeries_map_no_input_spec[this_str] = key

    def __init__(self, main_dict = None, shell: Shell = None, long_stf: Structure = None, ring_stf: Structure = None,
                 ring_frame: Structure = None):
        super(CylinderAndCurvedPlate, self).__init__()

        if all([main_dict is None, shell is None, long_stf is None, ring_stf is None]):
                self._sasd= None
                self._smsd= None
                self._tTsd= None
                self._tQsd= None
                self._psd = None
                self._shsd= None
                self._geometry= None
                self._mat_factor= None
                self._delta0 = None
                self._fab_method_ring_stf= None
                self._fab_method_ring_girder= None
                self._E = None
                self._v = None
                self._mat_yield= None
                self._length_between_girders= None
                self._panel_spacing= None
                self.__ring_stiffener_excluded= None
                self.__ring_frame_excluded= None
                self._end_cap_pressure_included= None
                self._uls_or_als= None
                self._cone_Nsd = None
                self._cone_M1sd = None
                self._cone_M2sd = None
                self._cone_Tsd = None
                self._cone_Q1sd = None
                self._cone_Q2sd = None

        else:
            self._sasd = main_dict['sasd'][0]
            self._smsd = main_dict['smsd'][0]
            self._tTsd = abs(main_dict['tTsd'][0])
            self._tQsd= main_dict['tQsd'][0]
            self._psd = main_dict['psd'][0]
            self._shsd = main_dict['shsd'][0]
            self._geometry = main_dict['geometry'][0]
            self._mat_factor = main_dict['material factor'][0]
            self._delta0 = main_dict['delta0'][0]
            self._fab_method_ring_stf = main_dict['fab method ring stf'][0]
            self._fab_method_ring_girder = main_dict['fab method ring girder'][0]
            self._E = main_dict['E-module'][0]
            self._v = main_dict['poisson'][0]
            self._mat_yield = main_dict['mat_yield'][0]
            self._length_between_girders = main_dict['length between girders'][0]
            self._panel_spacing = main_dict['panel spacing, s'][0]
            self.__ring_stiffener_excluded = main_dict['ring stf excluded'][0]
            self.__ring_frame_excluded = main_dict['ring frame excluded'][0]
            self._end_cap_pressure_included = main_dict['end cap pressure'][0]
            self._uls_or_als =  main_dict['ULS or ALS'][0]
            self._cone_Nsd = main_dict.get('cone Nsd', [None, 'kN'])[0]
            self._cone_M1sd = main_dict.get('cone M1sd', [None, 'kNm'])[0]
            self._cone_M2sd = main_dict.get('cone M2sd', [None, 'kNm'])[0]
            self._cone_Tsd = main_dict.get('cone Tsd', [None, 'kNm'])[0]
            self._cone_Q1sd = main_dict.get('cone Q1sd', [None, 'kN'])[0]
            self._cone_Q2sd = main_dict.get('cone Q2sd', [None, 'kN'])[0]

        self._Shell = shell
        self._LongStf = long_stf
        self._RingStf = ring_stf
        self._RingFrame = ring_frame



    def __str__(self):
        '''
        Returning all properties.
        '''

        long_string = 'N/A' if self._LongStf is None else self._LongStf.get_beam_string()
        ring_string = 'N/A' if self._RingStf is None else self._RingStf.get_beam_string()
        frame_string = 'N/A' if self._RingFrame is None else self._RingFrame.get_beam_string()
        s = max([self._Shell.dist_between_rings, 2*math.pi*self._Shell.radius])*1000 if self._LongStf == None else \
            self._LongStf.spacing

        return \
            str(
            '\n Cylinder radius:               ' + str(round(self._Shell.radius,3)) + ' meters' +
            '\n Cylinder thickness:            ' + str(self._Shell.thk*1000)+' mm'+
            '\n Distance between rings, l:     ' + str(self._Shell.dist_between_rings*1000)+' mm'+
            '\n Length of shell, L:            ' + str(self._Shell.length_of_shell*1000)+' mm'+
            '\n Total cylinder lenght:         ' + str(self._Shell.tot_cyl_length*1000)+' mm'+
            '\n Eff. Buckling length factor:   ' + str(self._Shell.k_factor)+
            '\n Material yield:                ' + str(self._mat_yield/1e6)+' MPa'+
            '\n Spacing/panel circ., s:        ' + str(s) + ' mm' +
            '\n Longitudinal stiffeners:       ' + long_string+
            '\n Ring stiffeners                ' + ring_string+
            '\n Ring frames/girders:           ' + frame_string+
            '\n Design axial stress/force:     ' + str(self._sasd/1e6)+' MPa'+
            '\n Design bending stress/moment:  ' + str(self._smsd/1e6)+' MPa'+
            '\n Design tosional stress/moment: ' + str(self._tTsd/1e6)+' MPa'+
            '\n Design shear stress/force:     ' + str(self._tQsd/1e6)+' MPa'+
            '\n Design lateral pressure        ' + str(self._psd/1e6)+' MPa'+
            '\n Additional hoop stress         ' + str(self._shsd/1e6)+' MPa')
    
    @property
    def uls_or_als(self):
        return self._uls_or_als
    @uls_or_als.setter
    def uls_or_als(self, val):
        self._uls_or_als = val
    @property
    def fab_method_ring_girder(self):
        return self._fab_method_ring_girder
    @fab_method_ring_girder.setter
    def fab_method_ring_girder(self, val):
        self._fab_method_ring_girder = val

    @property
    def fab_method_ring_stf(self):
        return self._fab_method_ring_stf
    @fab_method_ring_stf.setter
    def fab_method_ring_stf(self, val):
        self._fab_method_ring_stf = val
    @property
    def end_cap_pressure_included(self):
        return self._end_cap_pressure_included
    @end_cap_pressure_included.setter
    def end_cap_pressure_included(self, val):
        self._end_cap_pressure_included = val
    @property
    def delta0(self):
        return self._delta0
    @delta0.setter
    def delta0(self, val):
        self._delta0 = val
    @property
    def mat_factor(self):
        return self._mat_factor
    @mat_factor.setter
    def mat_factor(self, val):
        self._mat_factor = val
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
    def sasd(self):
        return self._sasd
    @sasd.setter
    def sasd(self, val):
        self._sasd = val
    @property
    def smsd(self):
        return self._smsd
    @smsd.setter
    def smsd(self, val):
        self._smsd = val
    @property
    def tTsd(self):
        return abs(self._tTsd)
    @tTsd.setter
    def tTsd(self, val):
        self._tTsd = abs(val)
    @property
    def tQsd(self):
        return self._tQsd
    @tQsd.setter
    def tQsd(self, val):
        self._tQsd = val
    @property
    def psd(self):
        return self._psd
    @psd.setter
    def psd(self, val):
        self._psd = val
    @property
    def shsd(self):
        return self._shsd
    @shsd.setter
    def shsd(self, val):
        self._shsd = val
    @property
    def panel_spacing(self):
        return self._panel_spacing
    @panel_spacing.setter
    def panel_spacing(self, val):
        self._panel_spacing = val
        
    @property
    def ShellObj(self):
        return self._Shell
    @ShellObj.setter
    def ShellObj(self, val):
        self._Shell = val
    @property
    def LongStfObj(self):
        return self._LongStf
    @LongStfObj.setter
    def LongStfObj(self, val):
        self._LongStf = val
    @property
    def RingStfObj(self):
        return self._RingStf
    @RingStfObj.setter
    def RingStfObj(self, val):
        self._RingStf = val
    @property
    def RingFrameObj(self):
        return self._RingFrame
    @RingFrameObj.setter
    def RingFrameObj(self, val):
        self._RingFrame = val

    @property
    def geometry(self):
        return self._geometry
    @geometry.setter
    def geometry(self, val):
        self._geometry = val
    @property
    def length_between_girders(self):
        return self._length_between_girders
    @length_between_girders.setter
    def length_between_girders(self, val):
        self._length_between_girders = val
    @property
    def _ring_stiffener_excluded(self):
        return self.__ring_stiffener_excluded
    @_ring_stiffener_excluded.setter
    def _ring_stiffener_excluded(self, val):
        self.__ring_stiffener_excluded = val
    @property
    def _ring_frame_excluded(self):
        return self.__ring_frame_excluded
    @_ring_frame_excluded.setter
    def _ring_frame_excluded(self, val):
        self.__ring_frame_excluded = val

    def get_utilization_factors(self, optimizing = False, empty_result_dict = False):
        '''
        If optimizing running time must be reduced.
        '''
        # Local buckling of stiffeners

        results = {'Unstiffened shell': None,
                   'Unstiffened conical shell': None,
                   'Unstiffened conical shell detailed': None,
                   'Longitudinal stiffened shell': None,
                   'Ring stiffened shell': None,
                   'Heavy ring frame': None,
                   'Column stability check': None,
                   'Column stability UF': None,
                   'Stiffener check': None,
                   'Stiffener check detailed': None,
                   'Weight': None}

        if empty_result_dict:
            return results
        if self._geometry == 9:
            conical = self.unstiffened_conical_shell()
            results['Unstiffened conical shell'] = conical['UF unstiffened conical shell']
            results['Unstiffened conical shell detailed'] = conical
            if optimizing:
                if results['Unstiffened conical shell'] > 1:
                    return False, 'UF unstiffened conical shell', results
                return True, 'Check OK', results
            return results
        data_shell_buckling = self.shell_buckling()
        unstiffend_shell, column_buckling_data = None, None
        # UF for unstiffened shell
        unstiffend_shell = self.unstiffened_shell(shell_data=data_shell_buckling)

        s = self._panel_spacing*1000 if self._LongStf is None else self._LongStf.spacing

        if any([self._geometry in [1, 5], s > self._Shell.dist_between_rings*1000]):
            uf_unstf_shell = unstiffend_shell['UF unstiffened circular cylinder']
            results['Unstiffened shell'] = uf_unstf_shell
        else:
            uf_unstf_shell = unstiffend_shell['UF unstiffened curved panel']
            results['Unstiffened shell'] = uf_unstf_shell

        if optimizing:
            if uf_unstf_shell > 1:
                return False, 'UF unstiffened', results

        # UF for longitudinal stiffened shell

        if self._geometry in [3,4,7,8]:
            if self._LongStf is not None:
                column_buckling_data= self.column_buckling(unstf_shell_data=unstiffend_shell,
                                                          shell_bukcling_data=data_shell_buckling)
                long_stf_shell = self.longitudinally_stiffened_shell(column_buckling_data=column_buckling_data,
                                                                     unstiffened_shell=unstiffend_shell)

                results['Column stability check'] = column_buckling_data['Column stability check']
                results['Column stability UF']  = column_buckling_data['Column stability UF']
                results['Need to check column buckling'] = column_buckling_data['Need to check column buckling']
                results['Stiffener check'] = column_buckling_data['stiffener check']
                results['Stiffener check detailed'] = column_buckling_data['stiffener check detailed']
                if self._geometry in [3,4,7,8] and long_stf_shell['fksd'] > 0:
                    results['Longitudinal stiffened shell'] = long_stf_shell['sjsd_used']/long_stf_shell['fksd']\
                        if self._geometry in [3,4,7,8] else 0

                if optimizing:
                    if not results['Column stability check']:
                        return False, 'Column stability', results
                    elif False in results['Stiffener check'].values():
                        return False, 'Stiffener check', results
                    elif results['Longitudinal stiffened shell'] > 1:
                        return False, 'UF longitudinal stiffeners', results

        if self._geometry in [5,6,7,8]:
            # UF for panel ring buckling
            ring_stf_shell = None
            if self._RingStf is not None:
                column_buckling_data = column_buckling_data if column_buckling_data is not None  \
                    else self.column_buckling( unstf_shell_data=unstiffend_shell,
                                               shell_bukcling_data=data_shell_buckling)
                ring_stf_shell = self.ring_stiffened_shell(data_shell_buckling=data_shell_buckling,
                                                           column_buckling_data=column_buckling_data)
                results['Column stability check'] = column_buckling_data['Column stability check']
                results['Column stability UF'] = column_buckling_data['Column stability UF']
                results['Need to check column buckling'] = column_buckling_data['Need to check column buckling']
                results['Stiffener check'] = column_buckling_data['stiffener check']
                results['Stiffener check detailed'] = column_buckling_data['stiffener check detailed']
                results['Ring stiffened shell'] = ring_stf_shell[0]



                if optimizing:
                    if not results['Column stability check']:
                        return False, 'Column stability', results
                    elif False in results['Stiffener check'].values():
                        return False, 'Stiffener check', results
                    elif results['Ring stiffened shell'] > 1:
                        return False, 'UF ring stiffeners', results

        # UF for ring frame
        if self._geometry in [5, 6, 7, 8]:
            if self._RingFrame is not None:
                column_buckling_data = column_buckling_data if column_buckling_data is not None  \
                    else self.column_buckling( unstf_shell_data=unstiffend_shell,
                                               shell_bukcling_data=data_shell_buckling)
                ring_stf_shell = ring_stf_shell if ring_stf_shell is not None else\
                    self.ring_stiffened_shell(data_shell_buckling=data_shell_buckling,
                                              column_buckling_data=column_buckling_data)
                results['Column stability check'] = column_buckling_data['Column stability check']
                results['Column stability UF'] = column_buckling_data['Column stability UF']
                results['Need to check column buckling'] = column_buckling_data['Need to check column buckling']
                results['Stiffener check'] = column_buckling_data['stiffener check']
                results['Stiffener check detailed'] = column_buckling_data['stiffener check detailed']
                results['Heavy ring frame'] = ring_stf_shell[1]

                if optimizing:
                    if not results['Column stability check']:
                        return False, 'Column stability', results
                    elif False in results['Stiffener check'].values():
                        return False, 'Stiffener check', results
                    elif results['Heavy ring frame'] > 1:
                        return False, 'UF ring frame', results

        if optimizing:
            return True, 'Check OK', results

        # print('Results for geometry', self._geometry)
        # print('UF',uf_unstf_shell, uf_long_stf, uf_ring_stf, uf_ring_frame)
        # print('Stiffeners', stiffener_check)

        return results

    def set_main_properties(self, main_dict):
        self._sasd = main_dict['sasd'][0]
        self._smsd = main_dict['smsd'][0]
        self._tTsd = abs(main_dict['tTsd'][0])
        self._tQsd= main_dict['tQsd'][0]
        self._psd = main_dict['psd'][0]
        self._shsd = main_dict['shsd'][0]
        self._geometry = main_dict['geometry'][0]
        self._mat_factor = main_dict['material factor'][0]
        self._delta0 = main_dict['delta0'][0]
        self._fab_method_ring_stf = main_dict['fab method ring stf'][0]
        self._fab_method_ring_girder = main_dict['fab method ring girder'][0]
        self._E = main_dict['E-module'][0]
        self._v = main_dict['poisson'][0]
        self._mat_yield = main_dict['mat_yield'][0]
        self._length_between_girders = main_dict['length between girders'][0]
        self._panel_spacing = main_dict['panel spacing, s'][0]
        self.__ring_stiffener_excluded = main_dict['ring stf excluded'][0]
        self.__ring_frame_excluded = main_dict['ring frame excluded'][0]
        self._end_cap_pressure_included = main_dict['end cap pressure'][0]
        self._uls_or_als =  main_dict['ULS or ALS'][0]
        self._cone_Nsd = main_dict.get('cone Nsd', [None, 'kN'])[0]
        self._cone_M1sd = main_dict.get('cone M1sd', [None, 'kNm'])[0]
        self._cone_M2sd = main_dict.get('cone M2sd', [None, 'kNm'])[0]
        self._cone_Tsd = main_dict.get('cone Tsd', [None, 'kNm'])[0]
        self._cone_Q1sd = main_dict.get('cone Q1sd', [None, 'kN'])[0]
        self._cone_Q2sd = main_dict.get('cone Q2sd', [None, 'kN'])[0]

    def shell_buckling(self):
        '''
        Main sheet to calculate cylinder buckling.
        '''
        stucture_objects = {'Unstiffened':self._Shell, 'Long Stiff.': self._LongStf, 'Ring Stiffeners': self._RingStf,
                            'Heavy ring Frame': self._RingFrame}
        stf_type = ['T', 'FB', 'T']
        assert self._Shell.dist_between_rings is not None, 'Input missing: self._Shell.dist_between_rings'
        assert self._Shell.radius is not None, 'Input missing: self._Shell.radius'
        assert self._Shell.thk is not None, 'Input missing: self._Shell.thk'

        l = self._Shell.dist_between_rings*1000
        r = self._Shell.radius*1000
        t = self._Shell.thk*1000

        parameters, cross_sec_data = list(), list()
        for idx, obj in stucture_objects.items():
            if obj is None:
                cross_sec_data.append([np.nan, np.nan, np.nan, np.nan, np.nan])
                if idx not in ['Unstiffened', 'Long Stiff.']:
                    parameters.append([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
                continue
            if idx != 'Unstiffened':
                hs = obj.hw/2 if obj.get_stiffener_type() == 'FB' else obj.hw + obj.tf/2
                It = obj.get_torsional_moment_venant()
                se = self._Shell.get_effective_width_shell_plate()

                Ipo = obj.get_polar_moment()
                Iz = obj.get_Iz_moment_of_inertia()

                Iy = obj.get_moment_of_intertia(efficent_se=se, tf1=self._Shell.thk)*1000**4

                cross_sec_data.append([hs, It, Iz, Ipo, Iy])

                A = obj.get_cross_section_area(include_plate=False)*math.pow(1000,2)
                beta = l/(1.56*math.sqrt(r*t))
                leo = (l/beta) *  ((math.cosh(2*beta)-math.cos(2*beta))/(math.sinh(2*beta)+math.sin(2*beta)))

                assert self._sasd is not None, 'Input missing: self._sasd'
                assert self._smsd is not None, 'Input missing: self._smsd'


                worst_axial_comb = min(self._sasd/1e6 - self._smsd/1e6, self._sasd/1e6 + self._smsd/1e6)
                sxsd_used = worst_axial_comb

                if idx == 'Long Stiff.':
                    zp = obj.get_cross_section_centroid_with_effective_plate(include_plate=False) * 1000
                    h_tot = obj.hw + obj.tf
                    zt = h_tot -zp
                else:
                    se = self._Shell.get_effective_width_shell_plate()
                    zp = obj.get_cross_section_centroid_with_effective_plate(se=se, tf1=self._Shell.thk) * 1000 # ch7.5.1 page 19
                    h_tot = self._Shell.thk*1000 + obj.hw + obj.tf
                    zt = h_tot -zp

            if idx not in ['Unstiffened', 'Long Stiff.']:  # Parameters
                alpha = A/(leo*t)
                zeta = max([0, 2*(math.sinh(beta)*math.cos(beta)+math.cosh(beta)*math.sin(beta))/
                            (math.sinh(2*beta)+math.sin(2*beta))])
                rf = r - t / 2 - (obj.hw + obj.tf)

                r0 = zt + rf
                parameters.append([alpha, beta, leo, zeta, rf, r0, zt])

        sxsd, shsd, shRsd, tsd = list(), list(), list(), list()

        for idx, obj in stucture_objects.items():
            if obj is None:
                shRsd.append(np.nan)
                continue
            assert self._psd is not None, 'Input missing: self._psd'

            if idx == 'Unstiffened':
                shsd.append((self._psd/1e6)*r/t + self._shsd/1e6)
                sxsd.append(self._sasd/1e6+self._smsd/1e6 if self._geometry in [2,6] else
                            min([self._sasd/1e6, self._sasd/1e6-self._smsd/1e6, self._sasd/1e6+self._smsd/1e6]))
                tsd.append(self._tTsd/1e6 + self._tQsd/1e6)
            elif idx == 'Long Stiff.':
                assert +self._shsd is not None, 'Input missing: self._shsd'
                assert +self._tTsd is not None, 'Input missing: self._tTsd'
                assert +self._tQsd is not None, 'Input missing: self._tQsd'

                if stucture_objects['Ring Stiffeners'] == None:
                    shsd.append(shsd[0]+self._shsd/1e6)
                else:
                    shsd_ring = ((self._psd/1e6)*r/t)-parameters[0][0]*parameters[0][3]/(parameters[0][0]+1)*\
                                ((self._psd/1e6)*r/t-0.3*sxsd[0])
                    shsd.append(shsd_ring + self._shsd/1e6)

                if self._geometry in [3,4,7,8]:
                    sxsd.append(sxsd_used)
                else:
                    sxsd.append(sxsd[0])

                tsd.append(self._tTsd/1e6 + self._tQsd/1e6)

            elif idx == 'Ring Stiffeners':
                rf = parameters[0][4]
                shsd_ring = ((self._psd / 1e6) * r / t) - parameters[0][0] * parameters[0][3] / (parameters[0][0] + 1) * \
                            ((self._psd / 1e6) * r / t - 0.3 * sxsd[0])
                shsd.append(np.nan if stucture_objects['Ring Stiffeners'] == None else shsd_ring)
                shRsd.append(((self._psd/1e6)*r/t-0.3*sxsd[0])*(1/(1+parameters[0][0]))*(r/rf))
                if self._geometry > 4:
                    sxsd.append(sxsd[0])
                    tsd.append(tsd[0])
                else:
                    sxsd.append(np.nan)
                    tsd.append(np.nan)

            else:
                rf = parameters[1][4]
                shsd.append(((self._psd/1e6)*r/t)-parameters[1][0]*parameters[1][3]/(parameters[1][0]+1)*
                            ((self._psd/1e6)*r/t-0.3*self._sasd/1e6))
                shRsd.append(((self._psd/1e6)*r/t-0.3*self._sasd/1e6)*(1/(1+parameters[1][0]))*(r/rf))
                if self._geometry > 4:
                    sxsd.append(sxsd[0])
                    tsd.append(tsd[0])
                else:
                    sxsd.append(np.nan)
                    tsd.append(np.nan)

        sxsd = np.array(sxsd)
        shsd = np.array(shsd)
        tsd = np.array(np.abs(tsd))
        sjsd = np.sqrt(sxsd**2 - sxsd*shsd + shsd**2+3*tsd**2)

        return {'sjsd': sjsd, 'parameters': parameters, 'cross section data': cross_sec_data,
                'shRsd': shRsd, 'shsd': shsd, 'sxsd': sxsd}

    def conical_stress_state(self, radius):
        assert self._Shell.cone_alpha is not None, 'Input missing: self._Shell.cone_alpha'
        assert self._Shell.thk is not None, 'Input missing: self._Shell.thk'
        cos_alpha = math.cos(math.radians(self._Shell.cone_alpha))
        te = self._Shell.thk * cos_alpha
        assert radius > 0, 'Input missing: conical shell radius'
        assert te > 0, 'Invalid conical shell equivalent thickness'

        psd = self._psd / 1e6
        nsd = getattr(self, '_cone_Nsd', 0) or 0
        m1sd = getattr(self, '_cone_M1sd', 0) or 0
        m2sd = getattr(self, '_cone_M2sd', 0) or 0
        tsd_force = getattr(self, '_cone_Tsd', 0) or 0
        q1sd = getattr(self, '_cone_Q1sd', 0) or 0
        q2sd = getattr(self, '_cone_Q2sd', 0) or 0

        pressure_pa = psd * 1e6
        axial = pressure_pa * radius / (2 * te) + nsd * 1000 / (2 * math.pi * radius * te)
        bending = math.sqrt(m1sd ** 2 + m2sd ** 2) * 1000 / (math.pi * radius ** 2 * te)
        hoop = pressure_pa * radius / te
        torsion = tsd_force * 1000 / (2 * math.pi * radius ** 2 * te)
        shear = math.sqrt(q1sd ** 2 + q2sd ** 2) * 1000 / (math.pi * radius * te)

        return {
            'radius': radius,
            'te': te,
            'sasd': axial,
            'smsd': bending,
            'shsd': hoop,
            'tTsd': torsion,
            'tQsd': shear,
        }

    def unstiffened_conical_shell(self):
        shell = self._Shell
        assert shell.cone_r1 is not None, 'Input missing: self._Shell.cone_r1'
        assert shell.cone_r2 is not None, 'Input missing: self._Shell.cone_r2'
        assert shell.cone_length is not None, 'Input missing: self._Shell.cone_length'

        original = {
            'geometry': self._geometry,
            'sasd': self._sasd,
            'smsd': self._smsd,
            'tTsd': self._tTsd,
            'tQsd': self._tQsd,
            'psd': self._psd,
            'shsd': self._shsd,
            'radius': shell.radius,
            'dist_between_rings': shell.dist_between_rings,
            'length_of_shell': shell.length_of_shell,
            'tot_cyl_length': shell.tot_cyl_length,
        }

        equivalent_radius = shell.cone_equivalent_radius()
        equivalent_length = shell.cone_equivalent_length()
        r_min = min(shell.cone_r1, shell.cone_r2)
        r_max = max(shell.cone_r1, shell.cone_r2)
        if r_min == r_max:
            radii = [r_min]
        else:
            radii = [r_min + (r_max - r_min) * idx / 100 for idx in range(101)]

        envelope = None
        try:
            for radius in radii:
                stress_state = self.conical_stress_state(radius)
                self._geometry = 1
                self._sasd = stress_state['sasd']
                self._smsd = stress_state['smsd']
                self._tTsd = abs(stress_state['tTsd'])
                self._tQsd = abs(stress_state['tQsd'])
                self._psd = 0
                self._shsd = stress_state['shsd']
                shell.radius = equivalent_radius
                shell.dist_between_rings = equivalent_length
                shell.length_of_shell = equivalent_length
                shell.tot_cyl_length = equivalent_length

                shell_data = self.shell_buckling()
                unstiffened = self.unstiffened_shell(shell_data=shell_data)
                uf = unstiffened['UF unstiffened circular cylinder']
                if envelope is None or uf > envelope['UF unstiffened conical shell']:
                    envelope = {
                        'UF unstiffened conical shell': uf,
                        'governing radius': radius,
                        'equivalent radius': equivalent_radius,
                        'equivalent length': equivalent_length,
                        'cone alpha': shell.cone_alpha,
                        'sasd': stress_state['sasd'],
                        'smsd': stress_state['smsd'],
                        'shsd': stress_state['shsd'],
                        'tTsd': stress_state['tTsd'],
                        'tQsd': stress_state['tQsd'],
                        'te': stress_state['te'],
                        'unstiffened shell data': unstiffened,
                    }
        finally:
            self._geometry = original['geometry']
            self._sasd = original['sasd']
            self._smsd = original['smsd']
            self._tTsd = original['tTsd']
            self._tQsd = original['tQsd']
            self._psd = original['psd']
            self._shsd = original['shsd']
            shell.radius = original['radius']
            shell.dist_between_rings = original['dist_between_rings']
            shell.length_of_shell = original['length_of_shell']
            shell.tot_cyl_length = original['tot_cyl_length']

        return envelope

    def unstiffened_shell(self, conical = False, shell_data = None):

        E = self._E/1e6
        t = self._Shell.thk*1000

        # get correct s

        s = min([self._Shell.dist_between_rings, 2*math.pi*self._Shell.radius])*1000 if self._LongStf == None else \
            self._LongStf.spacing
        v = self._v
        r = self._Shell.radius*1000
        l = self._Shell.dist_between_rings * 1000
        fy = self._mat_yield/1e6
        sasd = self._sasd/1e6
        smsd = self._smsd/1e6

        tsd = abs(self._tTsd/1e6+self._tQsd/1e6)

        psd = self._psd/1e6

        if self._RingStf is None:
            shsd = shell_data['shsd'][0]
        else:
            shsd = shell_data['shsd'][1]

        provide_data = dict()


        '''
        	Selections for: Type of Structure Geometry:
        1	Unstiffened shell (Force input)
        2	Unstiffened panel (Stress input)
        3	Longitudinal Stiffened shell  (Force input)
        4	Longitudinal Stiffened panel (Stress input)
        5	Ring Stiffened shell (Force input)
        6	Ring Stiffened panel (Stress input)
        7	Orthogonally Stiffened shell (Force input)
        8	Orthogonally Stiffened panel (Stress input)
        Selected:	
        3	Longitudinal Stiffened shell  (Force input)
        '''
        #   Pnt. 3.3 Unstifffed curved panel
        geometry = self._geometry

        if geometry in [2,6]:
            sxsd = sasd+smsd
        else:
            sxsd = min(sasd, sasd+smsd, sasd-smsd)

        if smsd < 0:
            smsd = abs(smsd)
            sm0sd = abs(smsd)

        else:
            if geometry in [2, 6]:
                smsd = 0
                sm0sd = 0
            else:
                smsd = smsd
                sm0sd = smsd

        sjsd = math.sqrt(math.pow(sxsd,2) - sxsd*shsd + math.pow(shsd,2) + 3 * math.pow(tsd, 2))  # (3.2.3)

        Zs = (math.pow(s, 2) / (r * t)) * math.sqrt(1 - math.pow(v, 2))  # The curvature parameter Zs (3.3.3)

        def table_3_1(chk):
            psi = {'Axial stress': 4, 'Shear stress': 5.34+4*math.pow(s/l, 2),
                   'Circumferential compression': math.pow(1+math.pow(s/l, 2), 2)}                      # ψ
            epsilon = {'Axial stress': 0.702*Zs, 'Shear stress': 0.856*math.sqrt(s/l)*math.pow(Zs, 3/4),
                   'Circumferential compression': 1.04*(s/l)*math.sqrt(Zs)}                             # ξ
            rho = {'Axial stress': 0.5*math.pow(1+(r/(150*t)), -0.5), 'Shear stress': 0.6,
                   'Circumferential compression': 0.6}
            return psi[chk], epsilon[chk], rho[chk]
        vals = list()

        for chk in ['Axial stress', 'Shear stress', 'Circumferential compression']:
            psi, epsilon, rho = table_3_1(chk=chk)
            C = psi * math.sqrt(1 + math.pow(rho * epsilon / psi, 2))  # (3.4.2) (3.6.4)
            fE = C*(math.pow(math.pi, 2)*E/(12*(1-math.pow(v,2)))) *math.pow(t/s,2)
            #print(chk, 'C', C, 'psi', psi,'epsilon', epsilon,'rho' ,rho, 'fE', fE)
            vals.append(fE)

        fEax, fEshear, fEcirc = vals
        sa0sd = -sxsd if sxsd < 0 else 0
        sh0sd = -shsd if shsd < 0 else 0 # Maximium allowable stress from iteration.

        if any([val == 0 for val in vals]):
            lambda_s_pow = 0
        else:
            lambda_s_pow = (fy/sjsd) * (sa0sd/fEax + sh0sd/fEcirc + tsd/fEshear)

        lambda_s = math.sqrt(lambda_s_pow)
        fks = fy/math.sqrt(1+math.pow(lambda_s,4 ))

        provide_data['fks - Unstifffed curved panel'] = fks
        if lambda_s < 0.5:
            gammaM = self._mat_factor
        else:
            if self._mat_factor == 1.1:
                if lambda_s > 1:
                    gammaM = 1.4
                else:
                    gammaM = 0.8+0.6*lambda_s
            elif self._mat_factor == 1.15:
                if lambda_s > 1:
                    gammaM = 1.45
                else:
                    gammaM = 0.85+0.6*lambda_s
            else:
                if lambda_s > 1:
                    gammaM = 1.45 * (self._mat_factor/1.15)
                else:
                    gammaM = 0.85+0.6*lambda_s * (self._mat_factor/1.15)
        if self._uls_or_als == 'ALS':
            gammaM = gammaM/self._mat_factor
        provide_data['gammaM Unstifffed panel'] = gammaM
        fksd = fks/gammaM
        provide_data['fksd - Unstifffed curved panel'] = fksd
        uf = sjsd/fksd

        provide_data['UF unstiffened curved panel'] = uf
        provide_data['gammaM curved panel'] = gammaM
        sjsd_max = math.sqrt(math.pow(sasd+smsd,2)-(sasd+smsd)*shsd+math.pow(shsd,2)+3*math.pow(tsd,2))

        uf_max =  self._mat_factor* sjsd_max/fy


        def iter_table_1():
            found, sasd_iter, count, this_val, logger  = False, 0.001 if uf > 1 else sasd, 0, 0, list()

            while not found:
                # Iteration
                sigmsd_iter = smsd if geometry in [2,6] else min([-smsd, smsd])
                siga0sd_iter = 0 if sasd_iter >= 0 else -sasd_iter  # (3.2.4)
                sigm0sd_iter = 0 if sigmsd_iter >= 0 else -sigmsd_iter  # (3.2.5)
                sigh0sd_iter = 0 if shsd>= 0 else -shsd  # (3.2.6)

                sjsd_iter = math.sqrt(math.pow(sasd_iter+sigmsd_iter, 2) - (sasd_iter+sigmsd_iter)*shsd + math.pow(shsd, 2)+
                                      3*math.pow(tsd, 2)) #(3.2.3)
                lambdas_iter = math.sqrt((fy / sjsd_iter) * ((siga0sd_iter+sigm0sd_iter)/fEax+ sigh0sd_iter/fEcirc+tsd/fEshear)) # (3.2.2)

                gammaM_iter = 1  # As taken in the DNVGL sheets
                fks_iter = fy / math.sqrt(1 + math.pow(lambdas_iter,4))
                fksd_iter = fks_iter / gammaM_iter
                #print('sjsd', sjsd_iter, 'fksd', fksd_iter, 'fks', fks, 'gammaM', gammaM_iter, 'lambdas_iter', lambdas_iter)
                this_val = sjsd_iter/fksd_iter
                logger.append(siga0sd_iter)
                if this_val > 1.0 or count == 1e6:
                    found = True
                count += 1
                if this_val >0.98:
                    sasd_iter -= 0.5
                elif this_val > 0.95:
                    sasd_iter -= 1
                elif this_val > 0.9:
                    sasd_iter -= 2
                elif this_val > 0.7:
                    sasd_iter -= 10
                else:
                    sasd_iter -= 20

                #print(sasd_iter, this_val)

            return 0 if len(logger) == 1 else max(max(logger[:-1]), 0)

        provide_data['max axial stress - 3.3 Unstifffed curved panel'] = iter_table_1()


        # Pnt. 3.4 Unstifffed circular cylinders
        Zl = (math.pow(l, 2)/(r*t)) * math.sqrt(1 - math.pow(v, 2)) #(3.4.3) (3.6.5)
        provide_data['Zl'] = Zl
        def table_3_2(chk):
            psi = {'Axial stress': 1, 'Bending': 1,
                   'Torsion and shear force': 5.34,
                   'Lateral pressure': 4, 'Hydrostatic pressure': 2}                      # ψ

            zeta= {'Axial stress': 0.702*Zl, 'Bending': 0.702*Zl,
                   'Torsion and shear force': 0.856* math.pow(Zl, 3/4),'Lateral pressure': 1.04*math.sqrt(Zl),
                   'Hydrostatic pressure': 1.04*math.sqrt(Zl)} # ξ

            rho = {'Axial stress': 0.5*math.pow(1+(r/(150*t)), -0.5), 'Bending': 0.5*math.pow(1+(r/(300*t)), -0.5),
                   'Torsion and shear force': 0.6,
                   'Lateral pressure': 0.6, 'Hydrostatic pressure': 0.6}
            return psi[chk], zeta[chk], rho[chk]

        vals = list()
        for chk in ['Axial stress', 'Bending', 'Torsion and shear force',
                    'Lateral pressure','Hydrostatic pressure']:
            psi, zeta, rho = table_3_2(chk=chk)
            C = psi * math.sqrt(1 + math.pow(rho * zeta / psi, 2))  # (3.4.2) (3.6.4)
            fE = C*math.pow(math.pi,2)*E / (12*(1-math.pow(v,2))) * math.pow(t/l,2)
            #print(chk, 'C', C, 'psi', psi,'epsilon', epsilon,'rho' ,rho, 'fE', fE)
            vals.append(fE)



        fEax, fEbend,  fEtors, fElat, fEhyd = vals

        provide_data['fEax - Unstifffed circular cylinders'] = fEax

        test1 = 3.85 * math.sqrt(r / t)
        test2 = 2.25 * math.sqrt(r / t)
        test_l_div_r = l/r
        provide_data['fEh - Unstifffed circular cylinders  - Psi=4'] = 0.25*E*math.pow(t/r,2) if test_l_div_r > test2 else fElat
        if l / r > test1:
            fEt_used = 0.25 * E * math.pow(t / r, 3 / 2)  # (3.4.4)
        else:
            fEt_used = fEtors

        if l / r > test2:
            fEh_used = 0.25 * E * math.pow(t / r, 2)
        else:
            fEh_used = fElat if self._end_cap_pressure_included == 'not included in axial stresses' else fEhyd

        sjsd = math.sqrt(math.pow(sxsd,2) - sxsd*shsd + math.pow(shsd,2) + 3 * math.pow(tsd, 2))  # (3.2.3)

        sa0sd = -sasd if sasd < 0 else 0
        sh0sd = -shsd if shsd < 0 else 0


        if any([fEax == 0, fEbend == 0, fEt_used == 0, fEh_used == 0, sjsd == 0]):
            lambda_s_pow = 0
        else:
            lambda_s_pow = (fy/sjsd) * (sa0sd/fEax + sm0sd/fEbend + sh0sd/fEh_used + tsd/fEt_used)

        lambda_s = math.sqrt(lambda_s_pow)
        fks = fy/math.sqrt(1+math.pow(lambda_s,4 ))

        provide_data['fks - Unstifffed circular cylinders'] = fks

        if lambda_s < 0.5:
            gammaM = self._mat_factor
        else:
            if self._mat_factor == 1.1:
                if lambda_s > 1:
                    gammaM = 1.4
                else:
                    gammaM = 0.8+0.6*lambda_s
            elif self._mat_factor == 1.15:
                if lambda_s > 1:
                    gammaM = 1.45
                else:
                    gammaM = 0.85+0.6*lambda_s
            else:
                if lambda_s > 1:
                    gammaM = 1.45 * (self._mat_factor/1.15)
                else:
                    gammaM = 0.85+0.6*lambda_s * (self._mat_factor/1.15)
        if self._uls_or_als == 'ALS':
            gammaM = gammaM/self._mat_factor

        fksd = fks/gammaM
        provide_data['fksd - Unstifffed circular cylinders'] = fksd
        uf = sjsd/fksd

        provide_data['UF unstiffened circular cylinder'] = uf
        provide_data['gammaM circular cylinder'] = gammaM
        #print('UF', uf, 'Unstifffed circular cylinders')
        def iter_table_2():
            found, sasd_iter, count, this_val, logger  = False, 0 if uf > 1 else sasd, 0, 0, list()
            while not found:
                # Iteration
                sigmsd_iter = smsd if geometry in [2, 6] else min([-smsd, smsd])
                siga0sd_iter = 0.00001 if sasd_iter >= 0 else -sasd_iter  # (3.2.4)
                sigm0sd_iter = 0.00001 if sigmsd_iter >= 0 else -sigmsd_iter  # (3.2.5)
                sigh0sd_iter = 0.00001 if shsd >= 0 else -shsd  # (3.2.6)
                sjsd_iter = math.sqrt(
                    math.pow(sasd_iter + sigmsd_iter, 2) - (sasd_iter + sigmsd_iter) * shsd + math.pow(shsd, 2) +
                    3 * math.pow(tsd, 2))  # (3.2.3)
                if sjsd_iter == 0:
                    sjsd_iter = 0.00001
                lambdas_iter = math.sqrt((fy/sjsd_iter) * (siga0sd_iter/fEax + sigm0sd_iter/fEbend +
                                                           sigh0sd_iter/fElat + tsd/fEtors))
                gammaM_iter = 1  # As taken in the DNVGL sheets
                fks_iter = fy / math.sqrt(1 + math.pow(lambdas_iter, 4))
                fksd_iter = fks_iter / gammaM_iter


                this_val = sjsd_iter / fksd_iter
                # print('sjsd', sjsd_iter, 'fksd', fksd_iter, 'fks', fks, 'gammaM', gammaM_iter, 'lambdas_iter',
                #       lambdas_iter, 'Count', count, 'UF', this_val)

                logger.append(siga0sd_iter)

                if this_val > 1.0 or count == 1e6:
                    found = True
                count += 1

                if this_val >0.98:
                    sasd_iter -= 0.5
                elif this_val > 0.95:
                    sasd_iter -= 1
                elif this_val > 0.9:
                    sasd_iter -= 2
                elif this_val > 0.7:
                    sasd_iter -= 10
                else:
                    sasd_iter -= 20
            return 0 if len(logger) == 1 else max(max(logger[:-1]),0)
        from_table = iter_table_2()

        provide_data['max axial stress - 3.4.2 Shell buckling'] = from_table

        provide_data['shsd'] = shsd
        return provide_data

    def ring_stiffened_shell(self, data_shell_buckling = None, column_buckling_data = None):

        E = self._E/1e6
        t = self._Shell.thk*1000
        s = min([self._Shell.dist_between_rings, 2*math.pi*self._Shell.radius])*1000 if self._LongStf == None else \
            self._LongStf.spacing

        r = self._Shell.radius*1000
        l = self._Shell.dist_between_rings * 1000
        fy = self._mat_yield/1e6

        L = self._Shell.tot_cyl_length*1000
        LH = L
        sasd = self._sasd/1e6
        smsd = self._smsd/1e6
        tsd = abs(self._tTsd/1e6 + self._tQsd/1e6)
        psd = self._psd/1e6

        data_shell_buckling = self.shell_buckling() if data_shell_buckling == None else data_shell_buckling

        #Pnt. 3.5:  Ring stiffened shell

        # Pnt. 3.5.2.1   Requirement for cross-sectional area:
        #Zl = self._Shell.get_Zl()

        Zl = math.pow(l, 2) * math.sqrt(1 - math.pow(self._v, 2)) / (r * t) if r * t > 0 else 0
        Areq = np.nan if Zl == 0 else (2/math.pow(Zl,2)+0.06)*l*t
        Areq = np.array([Areq, Areq])
        Astf = np.nan if self._RingStf is None else self._RingStf.get_cross_section_area(include_plate=False)*1000**2
        Aframe = np.nan if self._RingFrame is None else \
            self._RingFrame.get_cross_section_area(include_plate=False) * 1000 ** 2
        A = np.array([Astf, Aframe])

        uf_cross_section = Areq/A

        #Pnt. 3.5.2.3   Effective width calculation of shell plate
        lef = 1.56*math.sqrt(r*t)/(1+12*t/r)
        lef_used = np.array([min([lef, LH]), min([lef, LH])])

        #Pnt. 3.5.2.4   Required Ix for Shell subject to axial load
        A_long_stf = 0 if self._LongStf is None else self._LongStf.get_cross_section_area(include_plate=False)*1000**2
        alfaA = 0 if s*t <= 0 else A_long_stf/(s*t)


        r0 = np.array([data_shell_buckling['parameters'][0][5], data_shell_buckling['parameters'][1][5]])

        worst_ax_comp = min([sasd+smsd, sasd-smsd])

        Ixreq = np.array([abs(worst_ax_comp) * t * (1 + alfaA) * math.pow(r0[0], 4) / (500 * E * l),
                          abs(worst_ax_comp) * t * (1 + alfaA) * math.pow(r0[1], 4) / (500 * E * l)])

        #Pnt. 3.5.2.5   Required Ixh for shell subjected to torsion and/or shear:
        Ixhreq = np.array([math.pow(tsd / E, (8 / 5)) * math.pow(r0[0] / L, 1 / 5) * L * r0[0] * t * l,
                           math.pow(tsd / E, (8 / 5)) * math.pow(r0[1] / L, 1 / 5) * L * r0[1] * t * l])

        #Pnt. 3.5.2.6   Simplified calculation of Ih for shell subjected to external pressure
        zt = np.array([data_shell_buckling['parameters'][0][6],data_shell_buckling['parameters'][1][6]])
        rf = np.array([data_shell_buckling['parameters'][0][4], data_shell_buckling['parameters'][1][4]])

        delta0 = r*self._delta0

        fb_ring_req_val = np.array([0 if self._RingStf is None else 0.4*self._RingStf.tw*math.sqrt(E/fy),
                                    0 if self._RingFrame is None else 0.4*self._RingFrame.tw*math.sqrt(E/fy)])
        # if self._RingStf.get_stiffener_type() == 'FB':
        #     fb_ring_req = fb_ring_req_val[0] > self._RingStf.hw
        # else:
        #     fb_ring_req = np.NaN

        flanged_rf_req_h_val = np.array([0 if self._RingStf is None else 1.35*self._RingStf.tw*math.sqrt(E/fy),
                                         0 if self._RingFrame is None else 1.35*self._RingFrame.tw*math.sqrt(E/fy)])
        # if self._RingFrame.get_stiffener_type() != 'FB':
        #     flanged_rf_req_h = flanged_rf_req_h_val[1] > self._RingFrame.hw
        # else:
        #     flanged_rf_req_h = np.NaN

        flanged_rf_req_b_val = np.array([0 if self._RingStf is None else 7*self._RingStf.hw/math.sqrt(10+E*self._RingStf.hw/(fy*r)),
                                         0 if self._RingFrame is None else 7*self._RingFrame.hw/math.sqrt(10+E*self._RingFrame.hw/(fy*r))])
        # if self._RingFrame.get_stiffener_type() != 'FB':
        #     flanged_rf_req_b = flanged_rf_req_b_val[1] > self._RingFrame.b
        # else:
        #     flanged_rf_req_b = np.NaN

        if self._RingStf is not None:
            spf_stf = self._RingStf.hw/fb_ring_req_val[0] if self._RingStf.get_stiffener_type() == 'FB' \
                else max([flanged_rf_req_b_val[0]/self._RingStf.b, self._RingStf.hw/flanged_rf_req_h_val[0]])
        else:
            spf_stf = 0

        if self._RingFrame is not None:
            spf_frame = self._RingFrame.hw / fb_ring_req_val[1]if self._RingFrame.get_stiffener_type() == 'FB' \
                else max([flanged_rf_req_b_val[1] / self._RingFrame.b,self._RingFrame.hw / flanged_rf_req_h_val[1]])
        else:
            spf_frame = 0

        Stocky_profile_factor = np.array([spf_stf, spf_frame])

        fT = column_buckling_data['fT_dict']
        fT = np.array([fT['Ring Stiff.'] if Stocky_profile_factor[0] > 1 else fy,
                       fT['Ring Girder'] if Stocky_profile_factor[1] > 1 else fy])

        fr_used = np.array([fT[0] if self._fab_method_ring_stf == 1 else 0.9 * fT[0],
                            fT[1] if self._fab_method_ring_girder == 1 else 0.9 * fT[1]])
        shRsd = [abs(val) for val in data_shell_buckling['shRsd']]

        Ih = np.array([0 if E*r0[idx]*(fr_used[idx]/2-shRsd[idx]) == 0 else abs(psd)*r*math.pow(r0[idx],2)*l/(3*E)*
                                                                        (1.5+3*E*zt[idx]*delta0/(math.pow(r0[idx],2)
                                                                         *(fr_used[idx]/2-shRsd[idx])))
              for idx in [0,1]])

        # Pnt. 3.5.2.2     Moment of inertia:
        IR = [Ih[idx] + Ixhreq[idx] + Ixreq[idx] if all([psd <= 0, Ih[idx] > 0]) else Ixhreq[idx] + Ixreq[idx]
              for idx in [0,1]]
        Iy = [data_shell_buckling['cross section data'][idx+1][4] for idx in [0,1]]

        uf_moment_of_inertia = list()
        for idx in [0,1]:

            if Iy[idx] > 0:
                uf_moment_of_inertia.append(9.999 if fr_used[idx] < 2*shRsd[idx] else IR[idx]/Iy[idx])
            else:
                uf_moment_of_inertia.append(0)

        # Pnt. 3.5.2.7   Refined calculation of external pressure
        # parameters.append([alpha, beta, leo, zeta, rf, r0, zt])
        I = Iy
        Ihmax = [max(0, I[idx]-Ixhreq[idx]-Ixreq[idx]) for idx in [0,1]]
        leo = [data_shell_buckling['parameters'][idx][2] for idx in [0,1]]
        Ar = A
        ih2 = [0 if Ar[idx]+leo[idx]*t == 0 else Ihmax[idx]/(Ar[idx]+leo[idx]*t) for idx in [0,1]]
        alfa = [0 if l*t == 0 else 12*(1-math.pow(0.3,2))*Ihmax[idx]/(l*math.pow(t,3)) for idx in [0,1]]
        betta = [data_shell_buckling['parameters'][idx][0] for idx in [0,1]]
        ZL = [math.pow(L,2)/r/t*math.sqrt(1-math.pow(0.3,2)) for idx in [0,1]]

        C1 = [2*(1+alfa[idx])/(1+betta[idx])*(math.sqrt(1+0.27*ZL[idx]/math.sqrt(1+alfa[idx]))-alfa[idx]/(1+alfa[idx]))
              for idx in [0,1]]

        C2 = [2*math.sqrt(1+0.27*ZL[idx]) for idx in [0,1]]

        my = [0 if ih2[idx]*r*leo[idx]*C1[idx] == 0 else
              zt[idx]*delta0*rf[idx]*l/(ih2[idx]*r*leo[idx])*(1-C2[idx]/C1[idx])*1/(1-0.3/2) for idx in [0,1]]

        fE = np.array([C1[idx]*math.pow(math.pi,2)*E/(12*(1-math.pow(0.3,2)))*(math.pow(t/L,2)) if L > 0
                       else 0.1 for idx in [0,1]])

        fr = np.array(fT)
        lambda_2 = fr/fE
        lambda_ = np.sqrt(lambda_2)

        fk = [0 if lambda_2[idx] == 0 else fr[idx]*(1+my[idx]+lambda_2[idx]-math.sqrt(math.pow(1+my[idx]+lambda_2[idx],2)-
                                                                                 4*lambda_2[idx]))/(2*lambda_2[idx])
              for idx in [0,1]]
        gammaM = self._mat_factor # LRFD
        fkd = [fk[idx]/gammaM for idx in [0,1]]
        psd = np.array([0.75*fk[idx]*t*rf[idx]*(1+betta[idx])/(gammaM*math.pow(r,2)*(1-0.3/2)) for idx in [0,1]])

        uf_refined = abs((self._psd/1e6))/psd

        return np.max([uf_cross_section, uf_moment_of_inertia, uf_refined], axis=0)

    def longitudinally_stiffened_shell(self, column_buckling_data = None, unstiffened_shell = None):

        h = self._Shell.thk*1000 + self._LongStf.hw + self._LongStf.tf

        hw = self._LongStf.hw
        tw = self._LongStf.tw
        b = self._LongStf.b
        tf = self._LongStf.tf

        E = self._E/1e6
        t = self._Shell.thk*1000
        s = max([self._Shell.dist_between_rings, 2*math.pi*self._Shell.radius])*1000 if self._LongStf == None else \
            self._LongStf.spacing
        v = self._v
        r = self._Shell.radius*1000
        l = self._Shell.dist_between_rings * 1000
        fy = self._mat_yield/1e6

        L = self._Shell.tot_cyl_length*1000
        LH = L
        sasd = self._sasd/1e6
        smsd = self._smsd/1e6
        tsd = abs(self._tTsd/1e6 + self._tQsd/1e6)
        psd = self._psd/1e6
        shsd = unstiffened_shell['shsd']

        #print(h, hw, tw, b, tf, s, r, l, L, sasd, smsd, shsd)
        lightly_stf = s/t > math.sqrt(r/t)
        provide_data = dict()

        '''
        	Selections for: Type of Structure Geometry:
        1	Unstiffened shell (Force input)
        2	Unstiffened panel (Stress input)
        3	Longitudinal Stiffened shell  (Force input)
        4	Longitudinal Stiffened panel (Stress input)
        5	Ring Stiffened shell (Force input)
        6	Ring Stiffened panel (Stress input)
        7	Orthogonally Stiffened shell (Force input)
        8	Orthogonally Stiffened panel (Stress input)
        Selected:	
        3	Longitudinal Stiffened shell  (Force input)
        '''
        #   Pnt. 3.3 Unstifffed curved panel
        geometry = self._geometry
        data = unstiffened_shell if unstiffened_shell is not None else self.unstiffened_shell()

        if geometry == 1:
            fks = data['fks - Unstifffed circular cylinders']
        else:
            fks = data['fks - Unstifffed curved panel']

        sxSd  =min([sasd+smsd, sasd-smsd])


        sjsd  = math.sqrt(math.pow(sxSd,2) - sxSd*shsd + math.pow(shsd,2) + 3 * math.pow(tsd, 2))

        Se = (fks*abs(sxSd) / (sjsd*fy))*s

        # Moment of inertia
        As = A = hw*tw + b*tf  # checked

        num_stf = math.floor(2*math.pi*r/s)

        e= (hw*tw*(hw/2) + b*tf*(hw+tf/2)) / (hw*tw+b*tf)
        Istf = h*math.pow(tw,3)/12 + tf*math.pow(b, 3)/12

        dist_stf = r - t / 2 - e
        Istf_tot = 0
        angle = 0
        for stf_no in range(num_stf):
            Istf_tot += Istf + As*math.pow(dist_stf*math.cos(angle),2)
            angle += 2*math.pi/num_stf
        # Ishell = (math.pi/4) * ( math.pow(r+t/2,4) - math.pow(r-t/2,4))
        # Itot = Ishell + Istf_tot # Checked

        Iy = self._LongStf.get_moment_of_intertia(efficent_se=Se/1000, tf1=self._Shell.thk)*1000**4

        alpha = 12*(1-math.pow(v,2))*Iy/(s*math.pow(t,3))
        Zl = (math.pow(l, 2)/(r*t)) * math.sqrt(1-math.pow(v,2))

        #|1print('Zl', Zl, 'alpha', alpha, 'Isef', Iy, 'Se', Se, 'sjsd', sjsd, 'sxsd', sxSd, 'fks', fks, 'As', As)
        # Table 3-3


        def table_3_3(chk):
            psi = {'Axial stress': 0 if Se == 0 else (1+alpha) / (1+A/(Se*t)),
                   'Torsion and shear stress': 5.54+1.82*math.pow(l/s, 4/3) * math.pow(alpha, 1/3),
                   'Lateral Pressure': 2*(1+math.sqrt(1+alpha))}                      # ψ
            epsilon = {'Axial stress': 0.702*Zl,
                   'Torsion and shear stress': 0.856*math.pow(Zl, 3/4),
                   'Lateral Pressure': 1.04*math.sqrt(Zl)}                             # ξ
            rho = {'Axial stress': 0.5,
                   'Torsion and shear stress': 0.6,
                   'Lateral Pressure': 0.6}
            return psi[chk], epsilon[chk], rho[chk]

        vals = list()
        for chk in ['Axial stress', 'Torsion and shear stress','Lateral Pressure']:
            psi, epsilon, rho = table_3_3(chk=chk)

            C = 0 if psi == 0 else psi * math.sqrt(1 + math.pow(rho * epsilon / psi, 2))  # (3.4.2) (3.6.4)
            fE = C * ((math.pow(math.pi, 2) * E) / (12 * (1 - math.pow(v, 2)))) * math.pow(t / l,2)
            vals.append(fE)
            #print(chk, 'C', C, 'psi', psi,'epsilon', epsilon,'rho' ,rho, 'fE', fE)
        fEax, fEtors, fElat = vals

        #Torsional Buckling can be excluded as possible failure if:
        if self._LongStf._stiffener_type == 'FB':
            chk_fb = hw <= 0.4*tw*math.sqrt(E/fy)

        data_col_buc = column_buckling_data

        fy_used = fy if data_col_buc['lambda_T'] <= 0.6 else data_col_buc['fT']

        sasd = sasd*(A+s*t)/(A+Se*t) if A+Se*t>0 else 0
        smsd = smsd * (A + s * t) / (A + Se * t) if A + Se * t > 0 else 0

        sa0sd = -sasd if sasd < 0 else 0
        sm0sd = -smsd if smsd < 0 else 0
        sh0sd = -shsd if shsd < 0 else 0
        #print('fy_used', fy_used,'sasd', sasd,'shsd', shsd, 'tsd', tsd)
        sjsd_panels = math.sqrt(math.pow(sasd+smsd,2)-(sasd+smsd)*shsd + math.pow(shsd,2)+  3*math.pow(tsd,2))

        worst_axial_comb = min(sasd-smsd,sasd+smsd)
        sjsd_shells = math.sqrt(math.pow(worst_axial_comb,2)-worst_axial_comb*shsd +math.pow(shsd,2)+3*math.pow(tsd,2))
        sxsd_used = worst_axial_comb
        provide_data['sxsd_used'] = sxsd_used
        sjsd_used = sjsd_panels if self._geometry in [2,6] else sjsd_shells
        provide_data['sjsd_used'] = sjsd_used

        lambda_s2_panel = fy_used/sjsd_panels*((sa0sd+sm0sd)/fEax+sh0sd/fElat+tsd/fEtors) if\
            sjsd_panels*fEax*fEtors*fElat>0 else 0

        lambda_s2_shell = fy_used/sjsd_shells*(max(0,-worst_axial_comb)/fEax+sh0sd/fElat+tsd/fEtors) if\
            sjsd_shells*fEax*fEtors*fElat>0 else 0

        shell_type = 2 if self._geometry in [1,5] else 1
        lambda_s = math.sqrt(lambda_s2_panel) if shell_type == 1 else math.sqrt(lambda_s2_shell)

        fks = fy_used/math.sqrt(1+math.pow(lambda_s,4))

        #print('tsd',tsd, 'sasd', sasd, 'sjsd panels', sjsd_panels, 'fy_used', fy_used, 'lambda_T',data_col_buc['lambda_T'] )
        if lambda_s < 0.5:
            gammaM = self._mat_factor
        else:
            if self._mat_factor == 1.1:
                if lambda_s > 1:
                    gammaM = 1.4
                else:
                    gammaM = 0.8+0.6*lambda_s
            elif self._mat_factor == 1.15:
                if lambda_s > 1:
                    gammaM = 1.45
                else:
                    gammaM = 0.85+0.6*lambda_s
            else:
                if lambda_s > 1:
                    gammaM = 1.45 * (self._mat_factor/1.15)
                else:
                    gammaM = 0.85+0.6*lambda_s * (self._mat_factor/1.15)

        if self._uls_or_als == 'ALS':
            gammaM = gammaM/self._mat_factor

        # Design buckling strength:
        fksd = fks/gammaM
        provide_data['fksd'] = fksd
        # print(sjsd_panels, sjsd_shells, fksd)
        # print('fksd', fksd, 'fks', fks, 'gammaM', gammaM, 'lambda_s', lambda_s, 'lambda_s^2 panel',
        #       lambda_s2_panel, 'sjsd', sjsd_used, 'worst_axial_comb',worst_axial_comb, 'sm0sd',sm0sd)
        #print('  ')
        return provide_data

    @staticmethod
    def get_Itot(hw, tw, b, tf, r, s, t):

        h = t+hw+tf
        As = hw*tw + b*tf  # checked
        if As != 0:

            num_stf = math.floor(2*math.pi*r/s)
            e= (hw*tw*(hw/2) + b*tf*(hw+tf/2)) / (hw*tw+b*tf)
            Istf = h*math.pow(tw,3)/12 + tf*math.pow(b, 3)/12
            dist_stf = r - t / 2 - e
            Istf_tot = 0
            angle = 0
            for stf_no in range(num_stf):
                Istf_tot += Istf + As*math.pow(dist_stf*math.cos(angle),2)
                angle += 2*math.pi/num_stf
        else:
            Istf_tot = 0
        Ishell = (math.pi/4) * ( math.pow(r+t/2,4) - math.pow(r-t/2,4))
        Itot = Ishell + Istf_tot # Checked

        return Itot

    @staticmethod
    def _column_reduced_buckling_strength(fak, lambda_value):
        return (1 - 0.28 * math.pow(lambda_value, 2)) * fak if lambda_value <= 1.34 \
            else fak / math.pow(lambda_value, 2)

    def column_buckling(self,shell_bukcling_data = None, unstf_shell_data = None):

        geometry = self._geometry
        provide_data = dict()
        G = 80769.2

        if self._LongStf is None:
            h = self._Shell.thk*1000
        else:
            h = self._Shell.thk*1000 + self._LongStf.hw + self._LongStf.tf

        hw = 0 if self._LongStf is None else self._LongStf.hw
        tw = 0 if self._LongStf is None else self._LongStf.tw
        b = 0 if self._LongStf is None else self._LongStf.b
        tf = 0 if self._LongStf is None else self._LongStf.tf

        E = self._E/1e6
        t = self._Shell.thk*1000
        s = max([self._Shell.dist_between_rings, 2*math.pi*self._Shell.radius])*1000 if self._LongStf == None else \
            self._LongStf.spacing
        v = self._v
        r = self._Shell.radius*1000
        l = self._Shell.dist_between_rings * 1000
        fy = self._mat_yield/1e6

        L = self._Shell.tot_cyl_length*1000
        LH = L
        Lc = max([L, LH])

        sasd = self._sasd/1e6
        smsd = self._smsd/1e6
        tsd = abs(self._tTsd/1e6 + self._tQsd/1e6)
        psd = self._psd/1e6
        shsd = psd * r / t
        #print(t,h, hw, tw, b, tf, s, r, l, L, sasd, smsd, tsd, shsd)
        shell_buckling_data = self.shell_buckling() if\
            shell_bukcling_data is None else shell_bukcling_data
        data = self.unstiffened_shell(shell_data=shell_buckling_data) if unstf_shell_data is None else unstf_shell_data

        idx = 1
        param_map = {'Ring Stiff.': 0,'Ring Girder': 1}
        fT_dict = dict()
        for key, obj in {'Longitudinal stiff.': self._LongStf, 'Ring Stiff.': self._RingStf,
                         'Ring Girder': self._RingFrame}.items():
            if obj is None:
                idx += 1
                continue
            gammaM = data['gammaM circular cylinder'] if self._geometry > 2 else \
                data['gammaM curved panel']
            sjsd = shell_buckling_data['sjsd'][idx-1]

            this_s = 0 if self._LongStf is None else self._LongStf.spacing
            if any([self._geometry in [1, 5], this_s > (self._Shell.dist_between_rings * 1000)]):
                fksd = data['fksd - Unstifffed circular cylinders']
            else:
                fksd = data['fksd - Unstifffed curved panel']

            fks = fksd * gammaM
            eta = sjsd/fks
            hw = obj.hw
            tw = obj.tw
            b = obj.b
            tf = obj.tf

            if key == 'Longitudinal stiff.':
                s_or_leo = obj.spacing
                lT = l
            else:
                s_or_leo = shell_buckling_data['parameters'][param_map[key]][2]
                lT = math.pi*math.sqrt(r*hw)

            C = hw/s_or_leo*math.pow(t/tw,3)*math.sqrt(1-min([1,eta])) if s_or_leo*tw>0 else 0

            beta = (3*C+0.2)/(C+0.2)

            #parameters.append([alpha, beta, leo, zeta])
            hs, It, Iz, Ipo, Iy = shell_buckling_data['cross section data'][idx - 1]
            if obj.get_stiffener_type() == 'FB':
                Af = obj.tf * obj.b
                Aw = obj.hw * obj.tw
                fEt = beta * (Aw + math.pow(obj.tf / obj.tw, 2) * Af) / (Aw + 3 * Af) * G * math.pow(obj.tw / hw,
                                                                                                     2) + math.pow(
                    math.pi, 2) \
                      * E * Iz / ((Aw / 3 + Af) * math.pow(lT, 2))

            else:
                hs, It, Iz, Ipo, Iy = shell_buckling_data['cross section data'][idx-1]
                fEt = beta * G * It / Ipo + math.pow(math.pi, 2) * E * math.pow(hs, 2) * Iz / (Ipo * math.pow(lT, 2))
                #print(key, 'hs', hs, 'It', It, 'Iz', Iz, 'Ipo', Ipo, 'Iy', Iy)

            lambdaT = math.sqrt(fy/fEt)

            mu = 0.35*(lambdaT-0.6)
            fT = (1+mu+math.pow(lambdaT,2)-math.sqrt(math.pow(1+mu+math.pow(lambdaT,2),2)-4*math.pow(lambdaT,2)))\
                 /(2*math.pow(lambdaT,2))*fy if lambdaT > 0.6 else fy

            # General

            if key == 'Longitudinal stiff.':
                #print('Column buckling', 'fET', fEt, 'mu', mu, 'lambdaT', lambdaT, 'hs', hs, 'It', It, 'Iz', Iz ,'Ipo', Ipo)
                provide_data['lambda_T'] = lambdaT
                provide_data['fT'] = fT
            fT_dict[key] = fT
            idx += 1
            # if key == 'Ring Stiff.':
            #     print(hs, It, Iz, Ipo, Iy)
            #     print('hello')
        provide_data['fT_dict'] = fT_dict

        # Moment of inertia
        As = A = hw*tw + b*tf  # checked

        num_stf = math.floor(2*math.pi*r/s)

        Atot = As*num_stf + 2*math.pi*r*t

        e= (hw*tw*(hw/2) + b*tf*(hw+tf/2)) / (hw*tw+b*tf)
        Istf = h*math.pow(tw,3)/12 + tf*math.pow(b, 3)/12

        dist_stf = r - t / 2 - e
        Istf_tot = 0
        angle = 0
        for stf_no in range(num_stf):
            Istf_tot += Istf + As*math.pow(dist_stf*math.cos(angle),2)
            angle += 2*math.pi/num_stf

        Ishell = (math.pi/4) * ( math.pow(r+t/2,4) - math.pow(r-t/2,4))
        Itot = Ishell + Istf_tot # Checked

        k_factor = self._Shell.k_factor
        col_test =math.pow(k_factor*Lc/math.sqrt(Itot/Atot),2) >= 2.5*E/fy

        provide_data['Need to check column buckling'] = col_test
        # print("Column buckling should be assessed") if col_test else \
        #     print("Column buckling does not need to be checked")


        #Sec. 3.8.2   Column buckling strength:

        fEa = data['fEax - Unstifffed circular cylinders']
        #fEa = any([geometry in [1,5], s > l])
        fEh = data['fEh - Unstifffed circular cylinders  - Psi=4']

        #   Special case:  calculation of fak for unstiffened shell:

        #   General case:

        use_fac = 1 if geometry < 3 else 2

        if use_fac == 1:
            a = 1 + math.pow(fy, 2) / math.pow(fEa, 2)
            b = ((2 * math.pow(fy, 2) / (fEa * fEh)) - 1) * shsd
            c = math.pow(shsd, 2) + math.pow(fy, 2) * math.pow(shsd, 2) / math.pow(fEh, 2) - math.pow(fy, 2)
            fak = 0 if b == 0 else (b + math.sqrt(math.pow(b, 2) - 4 * a * c)) / (2 * a)
        elif any([geometry in [1,5], s > l]):
            fak = data['max axial stress - 3.4.2 Shell buckling']
        else:
            fak = data['max axial stress - 3.3 Unstifffed curved panel']

        i = math.sqrt(Itot/Atot)
        fE = 0.0001 if Lc*k_factor == 0 else E*math.pow(math.pi*i  / (Lc * k_factor), 2)

        Lambda_ = 0 if fE == 0 else math.sqrt(fak/fE)


        fkc = self._column_reduced_buckling_strength(fak, Lambda_)
        gammaM = data['gammaM curved panel'] #self._mat_factor  # Check

        fakd = fak/gammaM
        fkcd = fkc/gammaM

        sa0sd = -sasd if sasd<0 else 0

        # if fakd*fkcd > 0:
        #     stab_uf = sa0sd/fkcd + (abs(smsd) / (1-sa0sd/fE))/fakd
        #     stab_chk = stab_uf <= 1
        # else:
        #     stab_uf = 10
        #     stab_chk = True
        # print('use_fac', use_fac, 'geometry', geometry,  'Lambda_', Lambda_, 'math.sqrt(', fak, '/', fE)
        # print(sa0sd, '/', fkcd, '+', '(abs(', smsd, ') / (1 - ', sa0sd, ' / ', fE, ')) / ', fakd)
        if fakd*fkcd == 0:
            stab_uf = 0
            stab_chk = False
        else:
            stab_uf = sa0sd / fkcd + (abs(smsd) / (1 - sa0sd / fE)) / fakd
            stab_chk = stab_uf <= 1

        #print("Stability requirement satisfied") if stab_chk else print("Not acceptable")
        # Sec. 3.9   Torsional buckling:  moved to the top

        # Stiffener check

        stf_req_h = list()
        for idx, obj in enumerate([self._LongStf, self._RingStf, self._RingFrame]):
            if obj is None:
                stf_req_h.append(np.nan)
            else:
                stf_req_h.append(0.4*obj.tw*math.sqrt(E/fy) if obj.get_stiffener_type() == 'FB'
                                 else 1.35*obj.tw*math.sqrt(E/fy))

        stf_req_h = np.array(stf_req_h)

        stf_req_b = list()
        for idx, obj in enumerate([self._LongStf, self._RingStf, self._RingFrame]):
            if obj is None:
                stf_req_b.append(np.nan)
            else:
                stf_req_b.append(np.nan if obj.get_stiffener_type() == 'FB' else 0.4*obj.tf*math.sqrt(E/fy))

        bf = list()
        for idx, obj in enumerate([self._LongStf, self._RingStf, self._RingFrame]):
            if obj is None:
                bf.append(np.nan)
            elif obj.get_stiffener_type() == 'FB':
                bf.append(obj.b)
            elif obj.get_stiffener_type() == 'T':
                bf.append((obj.b-obj.tw)/2)
            else:
                bf.append(obj.b-obj.tw)
        bf = np.array(bf)

        hw_div_tw = list()
        for idx, obj in enumerate([self._RingStf, self._RingFrame]):
            if obj is None:
                hw_div_tw.append(np.nan)
            else:
                hw_div_tw.append(obj.hw/obj.tw)
        hw_div_tw = np.array(hw_div_tw)

        #parameters - [alpha, beta, leo, zeta, rf, r0, zt]

        req_hw_div_tw = list()
        for idx, obj in enumerate([self._RingStf, self._RingFrame]):
            if obj is None:
                req_hw_div_tw.append(np.nan)
            else:
                #print(shell_buckling_data['parameters'][idx][4],obj.tw, obj.hw, E, obj.hw, obj.b, obj.tf, fy)
                to_append = np.nan if obj.b*obj.tf == 0 else 2/3*math.sqrt(shell_buckling_data['parameters'][idx][4]
                                                                               *(obj.tw*obj.hw)*E/
                                                                               (obj.hw*obj.b*obj.tf*fy))
                req_hw_div_tw.append(to_append)
        req_hw_div_tw = np.array(req_hw_div_tw)

        ef_div_tw = list()
        for idx, obj in enumerate([self._RingStf, self._RingFrame]):
            if obj is None:
                ef_div_tw.append(np.nan)
            else:
                ef_div_tw.append(obj.get_flange_eccentricity())
        ef_div_tw = np.array(ef_div_tw)

        ef_div_tw_req = list()
        for idx, obj in enumerate([self._RingStf, self._RingFrame]):
            if obj is None:
                ef_div_tw_req.append(np.nan)
            else:
                ef_div_tw_req.append(np.nan if obj.b*obj.tf == 0 else
                             1/3*shell_buckling_data['parameters'][idx][4]/obj.hw*obj.hw*obj.tw/(obj.b*obj.tf))
        ef_div_tw_req = np.array(ef_div_tw_req)

        #
        # print(stf_req_h , '>', np.array([np.nan if self._LongStf is None else self._LongStf.hw,
        #                                  np.nan if self._RingStf is None else self._RingStf.hw,
        #                                  np.nan if self._RingFrame is None else self._RingFrame.hw]))
        # print(stf_req_b , '>', bf)
        # print(hw_div_tw , '<', req_hw_div_tw)
        # print(ef_div_tw , '<', ef_div_tw_req)

        chk1 = stf_req_h>np.array([np.nan if self._LongStf is None else self._LongStf.hw,
                                  np.nan if self._RingStf is None else self._RingStf.hw,
                                  np.nan if self._RingFrame is None else self._RingFrame.hw])
        chk1 = [np.nan if np.isnan(val) else chk1[idx] for idx, val in enumerate(stf_req_h)]

        chk2 = stf_req_b > bf
        chk2 = [np.nan if np.isnan(val) else chk2[idx] for idx, val in enumerate(stf_req_b)]

        chk3= hw_div_tw < req_hw_div_tw
        chk3 = [np.nan if np.isnan(val) else chk3[idx] for idx, val in enumerate(req_hw_div_tw)]

        chk4 = ef_div_tw < ef_div_tw_req
        chk4 = [np.nan if np.isnan(val) else chk4[idx] for idx, val in enumerate(ef_div_tw_req)]

        provide_data['stiffener check'] = {'longitudinal':all([chk1[0], chk2[0]]),
                                           'ring stiffener': None if self._RingStf is None else
                                           all([chk1[1],chk2[1],chk3[0],chk4[0]]),
                                           'ring frame': None if self._RingFrame is None else
                                           True}
                                           #all([chk1[2],chk2[2],chk3[1],chk4[1]])} SKIP check for girders
        provide_data['stiffener check detailed'] = {'longitudinal':'Web height < ' + str(round(stf_req_h[0],1)) if not chk1[0]
        else '' + ' ' + 'flange width < ' +str(round(stf_req_b[0],1)) if not chk2[0] else ' ',
                                                   'ring stiffener': None if self._RingStf is None
                                                   else 'Web height < ' + str(round(stf_req_h[1],1)) if not chk1[1]
                                                   else '' + ' ' + 'flange width < ' +str(round(stf_req_b[1],1)) if not chk2[1]
                                                   else ' '  + ' ' + 'hw/tw >= ' + str(round(req_hw_div_tw[0],1))
                                                   if not chk3[0]
                                                   else ''+ ' ' + 'ef/tw >= ' + str(round(ef_div_tw_req[0],1))
                                                   if not chk4[0]
                                                   else '',
                                                   'ring frame': None if self._RingFrame is None
                                                   else  'Web height < ' + str(round(stf_req_h[2],1)) if not chk1[2]
                                                   else '' + ' ' + 'flange width < ' +str(round(stf_req_b[2],1)) if not chk2[2]
                                                   else ' '  + ' ' + 'hw/tw >= ' + str(round(req_hw_div_tw[1],1))
                                                   if not chk3[1]
                                                   else ''+ ' ' + 'ef/tw >= ' + str(round(ef_div_tw_req[1],1))
                                                   if not chk4[1]
                                                   else ''}

        provide_data['Column stability check'] = stab_chk
        provide_data['Column stability UF'] = stab_uf

        return provide_data

    def get_all_properties(self):
        all_data = {'Main class': self.get_main_properties(),
                    'Shell': self._Shell.get_main_properties(),
                    'Long. stf.': None if self._LongStf is None else self._LongStf.get_structure_prop(),
                    'Ring stf.': None if self._RingStf is None else self.RingStfObj.get_structure_prop(),
                    'Ring frame': None if self._RingFrame is None else self._RingFrame.get_structure_prop()}
        return all_data

    def set_all_properties(self, all_prop_dict): # TODO ensure that this is set when optimizing and saving.
        all_data = {'Main class': self.set_main_properties(all_prop_dict['Main class']),
                    'Shell': self._Shell.set_main_properties(all_prop_dict['Shell']),
                    'Long. stf.': None if self._LongStf is None else
                    self._LongStf.set_main_properties(all_prop_dict['Long. stf.']),
                    'Ring stf.': None if self._RingStf is None else
                    self.RingStfObj.set_main_properties(all_prop_dict['Ring stf.']),
                    'Ring frame': None if self._RingFrame is None else
                    self._RingFrame.set_main_properties(all_prop_dict['Ring frame'])}
        return all_data

    def get_main_properties(self):
        main_dict = {'sasd': [self._sasd, 'Pa'],
                     'smsd': [self._smsd, 'Pa'],
                     'tTsd': [abs(self._tTsd), 'Pa'],
                     'tQsd': [self._tQsd, 'Pa'],
                     'psd': [self._psd, 'Pa'],
                     'shsd': [self._shsd, 'Pa'],
                     'geometry': [self._geometry, ''],
                     'material factor': [self._mat_factor, ''],
                     'delta0': [self._delta0, ''],
                     'fab method ring stf': [self._fab_method_ring_stf, '-'],
                     'fab method ring girder': [self._fab_method_ring_girder, '-'],
                     'E-module': [self._E, 'Pa'],
                     'poisson': [self._v, '-'],
                     'mat_yield': [self._mat_yield, 'Pa'],
                     'length between girders': [self._length_between_girders, 'm'],
                     'panel spacing, s':  [self._panel_spacing, 'm'],
                     'ring stf excluded': [self.__ring_stiffener_excluded, ''],
                     'ring frame excluded': [self.__ring_frame_excluded, ''],
                     'end cap pressure': [self._end_cap_pressure_included, ''],
                     'ULS or ALS':[self._uls_or_als, ''],
                     'cone Nsd': [self._cone_Nsd, 'kN'],
                     'cone M1sd': [self._cone_M1sd, 'kNm'],
                     'cone M2sd': [self._cone_M2sd, 'kNm'],
                     'cone Tsd': [self._cone_Tsd, 'kNm'],
                     'cone Q1sd': [self._cone_Q1sd, 'kN'],
                     'cone Q2sd': [self._cone_Q2sd, 'kN']}

        return main_dict
        
    def set_stresses_and_pressure(self, val):
        self._sasd = val['sasd']
        self._smsd = val['smsd']
        self._tTsd = abs(val['tTsd'])
        self._tQsd= val['tQsd']
        self._psd = val['psd']
        self._shsd = val['shsd']
        
    def get_x_opt(self):
        '''
        shell       (0.02, 2.5, 5, 5, 10, nan, nan, nan),
        long        (0.875, nan, 0.3, 0.01, 0.1, 0.01, nan, stiffener_type)),
        ring        (nan, nan, 0.3, 0.01, 0.1, 0.01, nan, stiffener_type)),
        ring        (nan, nan, 0.7, 0.02, 0.2, 0.02, nan, stiffener_type))] 
        
        (self._spacing, self._plate_th, self._web_height, self._web_th, self._flange_width,
                self._flange_th, self._span, self._girder_lg, self._stiffener_type)
        '''
        if self._geometry == 9:
            shell = [
                self._Shell.thk,
                self._Shell.cone_equivalent_radius(),
                self._Shell.cone_equivalent_length(),
                self._Shell.cone_equivalent_length(),
                self._Shell.cone_equivalent_length(),
                self._Shell.cone_r1,
                self._Shell.cone_r2,
                self._Shell.cone_length,
            ]
        else:
            shell = [self._Shell.thk, self._Shell.radius, self._Shell.dist_between_rings, self._Shell.length_of_shell,
                     self._Shell.tot_cyl_length, np.nan, np.nan, np.nan]
        if self._LongStf is not None:
            long = [self._LongStf.spacing/1000, np.nan, self._LongStf.hw/1000, self._LongStf.tw/1000, self._LongStf.b/1000, 
                    self._LongStf.tf/1000, np.nan, self._LongStf.stiffener_type]
        else:
            long = [0 for dummy in range(8)]
        
        if self._RingStf is not None:
            ring_stf = [self._RingStf.s/1000, np.nan, self._RingStf.hw/1000, self._RingStf.tw/1000, self._RingStf.b/1000, 
                    self._RingStf.tf/1000, np.nan, self._RingStf.stiffener_type]
        else:
            ring_stf = [0 for dummy in range(8)]
        
        if self._RingFrame is not None:
            ring_fr = [self._RingFrame.s/1000, np.nan, self._RingFrame.hw/1000, self._RingFrame.tw/1000, self._RingFrame.b/1000, 
                    self._RingFrame.tf/1000, np.nan, self._RingFrame.stiffener_type]
        else:
            ring_fr = [0 for dummy in range(8)]

        return [shell, long, ring_stf, ring_fr]

