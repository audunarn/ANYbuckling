'''
Example input dictionaries for ANYbuckling tests.

Extracted from ANYstructure anystruct/example_data.py.
'''

structure_types = {'vertical': ['BBS', 'SIDE_SHELL', 'SSS'],
                         'horizontal': ['BOTTOM', 'BBT', 'HOPPER', 'MD'],
                         'non-wt': ['FRAME', 'GENERAL_INTERNAL_NONWT'],
                         'internals': ['INNER_SIDE', 'FRAME_WT', 'GENERAL_INTERNAL_WT',
                                       'INTERNAL_ZERO_STRESS_WT', 'INTERNAL_LOW_STRESS_WT']}

obj_dict = {'mat_yield': [355e6, 'Pa'], 'mat_factor': [1.10, ''],'span': [3.3, 'm'], 'spacing': [0.68, 'm'],
            'plate_thk': [0.025, 'm'],
            'stf_web_height': [0.250297358, 'm'], 'stf_web_thk': [0.012, 'm'], 'stf_flange_width': [0.052, 'm'],
            'stf_flange_thk': [0.029702642, 'm'], 'structure_type': ['BOTTOM', ''], 'plate_kpp': [1, ''],
            'stf_kps': [1, ''], 'stf_km1': [12, ''], 'stf_km2': [24, ''], 'stf_km3': [12, ''],
            'sigma_y1': [100, 'MPa'], 'sigma_y2': [100, 'MPa'], 'sigma_x2': [102.7, 'MPa'], 'sigma_x1': [102.7, 'MPa'],
            'tau_xy': [5, 'MPa'],
            'stf_type': ['T', ''], 'structure_types': [structure_types, ''], 'zstar_optimization': [True, ''],
            'puls buckling method':[1,''], 'puls boundary':['Int',''], 'puls stiffener end':['C',''],
            'puls sp or up':['SP',''], 'puls up boundary' :['SSSS',''], 'panel or shell': ['panel', ''],
            'pressure side': ['both sides', ''], 'girder_lg': [5, 'm']}

obj_dict_cyl_ring2 = {'mat_yield': [355e6, 'Pa'], 'mat_factor': [1.15, ''],'span': [5, 'm'], 'spacing': [0.7, 'm'],
                    'plate_thk': [0.020, 'm'],
                    'stf_web_height': [0.3, 'm'], 'stf_web_thk': [0.012, 'm'], 'stf_flange_width': [0.12, 'm'],
                    'stf_flange_thk': [0.02, 'm'], 'structure_type': ['BOTTOM', ''], 'plate_kpp': [1, ''],
                    'stf_kps': [1, ''], 'stf_km1': [12, ''], 'stf_km2': [24, ''], 'stf_km3': [12, ''],
                    'sigma_y1': [80, 'MPa'], 'sigma_y2': [80, 'MPa'], 'sigma_x2': [80, 'MPa'], 'sigma_x1': [80, 'MPa'], 'tau_xy': [5, 'MPa'],
                    'stf_type': ['T', ''], 'structure_types': [structure_types, ''], 'zstar_optimization': [True, ''],
                    'puls buckling method':[2,''], 'puls boundary':['Int',''], 'puls stiffener end':['C',''],
                    'puls sp or up':['SP',''], 'puls up boundary' :['SSSS',''] , 'panel or shell': ['shell', ''] }

obj_dict_heavy = {'mat_yield': [355e6, 'Pa'], 'mat_factor': [1.15, ''],'span': [3700, 'm'], 'spacing': [0.75, 'm'],
            'plate_thk': [0.018, 'm'],
            'stf_web_height': [0.500, 'm'], 'stf_web_thk': [0.0120, 'm'], 'stf_flange_width': [0.150, 'm'],
            'stf_flange_thk': [0.02, 'm'], 'structure_type': ['BOTTOM', ''], 'plate_kpp': [1, ''],
            'stf_kps': [1, ''], 'stf_km1': [12, ''], 'stf_km2': [24, ''], 'stf_km3': [12, ''],
            'sigma_y1': [80, 'MPa'], 'sigma_y2': [80, 'MPa'], 'sigma_x2': [80, 'MPa'], 'sigma_x1': [80, 'MPa'], 'tau_xy': [5, 'MPa'],
            'stf_type': ['T', ''], 'structure_types': [structure_types, ''], 'zstar_optimization': [True, ''],
                  'puls buckling method':[2,''], 'puls boundary':['Int',''], 'puls stiffener end':['C',''],
            'puls sp or up':['SP',''], 'puls up boundary' :['SSSS',''], 'panel or shell': ['panel', ''],
            'pressure side': ['both sides', '']}

obj_dict2 = {'mat_yield': [355e6, 'Pa'], 'mat_factor': [1.15, ''],'span': [4, 'm'], 'spacing': [0.7, 'm'],
            'plate_thk': [0.018, 'm'],
            'stf_web_height': [0.36, 'm'], 'stf_web_thk': [0.012, 'm'], 'stf_flange_width': [0.15, 'm'],
            'stf_flange_thk': [0.02, 'm'], 'structure_type': ['BOTTOM', ''], 'plate_kpp': [1, ''],
            'stf_kps': [1, ''], 'stf_km1': [12, ''], 'stf_km2': [24, ''], 'stf_km3': [12, ''],
            'sigma_y1': [100, 'MPa'], 'sigma_y2': [100, 'MPa'], 'sigma_x2': [80, 'MPa'], 'sigma_x1': [50, 'MPa'], 'tau_xy': [5, 'MPa'],
            'stf_type': ['T', ''], 'structure_types': [structure_types, ''], 'zstar_optimization': [True, ''],
             'puls buckling method':[2,''], 'puls boundary':['Int',''], 'puls stiffener end':['C',''],
            'puls sp or up':['SP',''], 'puls up boundary' :['SSSS',''] }

obj_dict_L = {'mat_yield': [355e6, 'Pa'], 'mat_factor': [1.15, ''], 'span': [3.6, 'm'], 'spacing': [0.82, 'm'],
            'plate_thk': [0.018, 'm'],
            'stf_web_height': [0.4, 'm'], 'stf_web_thk': [0.014, 'm'], 'stf_flange_width': [0.072, 'm'],
            'stf_flange_thk': [0.0439, 'm'], 'structure_type': ['BOTTOM', ''], 'plate_kpp': [0.5, ''],
            'stf_kps': [1, ''], 'stf_km1': [12, ''], 'stf_km2': [24, ''], 'stf_km3': [12, ''],
            'sigma_y1': [102, 'MPa'], 'sigma_y2': [106.9, 'MPa'], 'sigma_x2': [66.8, 'MPa'], 'sigma_x1': [66.8, 'MPa'], 'tau_xy': [20, 'MPa'],
            'stf_type': ['L', ''], 'structure_types': [structure_types, ''], 'zstar_optimization': [True, ''],
              'puls buckling method':[2,''], 'puls boundary':['Int',''], 'puls stiffener end':['C',''],
            'puls sp or up':['SP',''], 'puls up boundary' :['SSSS',''] , 'panel or shell': ['panel', ''], 'pressure side': ['both sides', ''] }

shell_dict = {'plate_thk': [20 / 1000, 'm'],
              'radius': [5000 / 1000, 'm'],
              'distance between rings, l': [700 / 1000, 'm'],
              'length of shell, L': [5000 / 1000, 'm'],
              'tot cyl length, Lc': [5000 / 1000, 'm'],
              'eff. buckling lenght factor': [1, ''],
              'mat_yield': [355 * 1e6, 'Pa'],
              }

shell_main_dict2 = {'sasd': [79.58 * 1e6, 'Pa'],
             'smsd': [31.89* 1e6, 'Pa'],
             'tTsd': [12.73* 1e6, 'Pa'],
             'tQsd': [4.77* 1e6, 'Pa'],
             'psd': [-0.2* 1e6, 'Pa'],
             'shsd': [0, 'Pa'],
             'geometry': [5, '-'],
             'material factor': [1.15, ''],
             'delta0': [0.005, ''],
             'fab method ring stf': [1, ''],
             'fab method ring girder': [1, ''],
             'E-module': [2.1e11, 'Pa'],
             'poisson': [0.3, '-'],
             'mat_yield': [355 * 1e6, 'Pa'],
                    'length between girders': [None, 'm'],
                    'panel spacing, s': [0.7, 'm'],
                    'ring stf excluded': [False, ''],
                    'ring frame excluded': [True, ''],
                    'end cap pressure': ['not included in axial stresses', ''],
                   'ULS or ALS': ['ULS', '']}

prescriptive_main_dict = dict()
prescriptive_main_dict['minimum pressure in adjacent spans']  = [None, '']
prescriptive_main_dict['material yield']  = [355e6, 'Pa']
prescriptive_main_dict['load factor on stresses']  = [1, '']
prescriptive_main_dict['load factor on pressure']  = [1, '']
prescriptive_main_dict['buckling method']  = ['ultimate', '']
prescriptive_main_dict['stiffener end support']  = ['Continuous', '']  # 'Continuous'
prescriptive_main_dict['girder end support']  = ['Continuous', '']  # 'Continuous'
prescriptive_main_dict['tension field']  = ['not allowed', '']  # 'not allowed'
prescriptive_main_dict['plate effective agains sigy']   = [True, ''] # True
prescriptive_main_dict['buckling length factor stf']  = [None, '']
prescriptive_main_dict['buckling length factor girder']  = [None, '']
prescriptive_main_dict['km3']  = [12, '']  # 12
prescriptive_main_dict['km2']  = [24, '']  # 24
prescriptive_main_dict['girder distance between lateral support']  = [None, '']
prescriptive_main_dict['stiffener distance between lateral support']  = [None, '']
prescriptive_main_dict['kgirder']  = [None, '']
prescriptive_main_dict['panel length, Lp']  = [None, '']
prescriptive_main_dict['pressure side']  = ['both sides', '']# either 'stiffener', 'plate', 'both'
prescriptive_main_dict['fabrication method stiffener'] =  ['welded', '']
prescriptive_main_dict['fabrication method girder'] =   ['welded', '']
prescriptive_main_dict['calculation domain'] = ['Flat plate, stiffened', '']
