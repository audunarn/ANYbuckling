'''
Stiffener section library utilities.

Transferred from ANYstructure (anystruct/helper.py). Reads stiffener
section definitions from xml, json or csv files. The package ships a
bulb/anglebar/T-bar/flatbar section table used for bulb-profile
conversion (see :func:`bulb_section_file`).
'''
import copy
import csv
import json
from pathlib import Path

_DATA_DIR = Path(__file__).parent / 'data'


def bulb_section_file() -> str:
    ''' Return the path to the packaged bulb/anglebar/T-bar/flatbar section table. '''
    return str(_DATA_DIR / 'bulb_anglebar_tbar_flatbar.csv')


def helper_read_section_file(files, obj = None, to_json = False, to_csv = None):
    ''' Read a section definition file (xml, json or csv).

    Returns a list of section property dicts, or - if ``obj`` is given -
    copies of the object(s) with the section properties applied.
    '''
    from xml.dom import minidom
    to_return_final, to_return, return_csv = list(),  dict(), list()
    if type(files) != list:
        files = [files,]
    for file in files:
        if file.endswith('xml'):
            xmldoc = minidom.parse(file)
            sectionlist = xmldoc.getElementsByTagName('section')
            sec_types = ('unsymmetrical_i_section', 'l_section', 'bar_section')

            for idx, sec_type in enumerate(sec_types):
                sec_type_get = xmldoc.getElementsByTagName(sec_type)
                if sec_types == []:
                    continue
                for item, itemdata in zip(sectionlist, sec_type_get):
                    if sec_type == sec_types[0]:
                        stf_web_h, stf_web_thk = 'h', 'tw'
                        stf_flange_width, stf_flange_thk  = 'bfbot', 'tfbot'
                        stiffener_type = 'T'
                        mult = 1/1000
                    elif sec_type == sec_types[1]:
                        stf_web_h, stf_web_thk = 'h', 'tw'
                        stf_flange_width, stf_flange_thk  = 'b', 'tf'
                        stiffener_type = 'L'
                        mult = 1/1000
                    elif sec_type == sec_types[2]:
                        stf_web_h, stf_web_thk = 'h', 'b'
                        stf_flange_width, stf_flange_thk  = None, None
                        stiffener_type = 'FB'
                        mult = 1 / 1000
                    section_name = item.getAttribute('name')
                    to_return[section_name] = {'stf_web_height': [float(itemdata.getAttribute(stf_web_h)) *mult, 'm'],
                                               'stf_web_thk': [float(itemdata.getAttribute(stf_web_thk)) *mult,'m'],
                                               'stf_flange_width': [0 if stf_flange_width is None else
                                               float(itemdata.getAttribute(stf_flange_width)) *mult,'m'],
                                               'stf_flange_thk': [0 if stf_flange_thk is None else
                                               float(itemdata.getAttribute(stf_flange_thk)) *mult, 'm'],
                                               'stf_type': [stiffener_type, '']}

                    return_csv.append([to_return[section_name][var][0] for var in ['stf_web_height', 'stf_web_thk',
                                                                                   'stf_flange_width', 'stf_flange_thk',
                                                                                   'stf_type']])

        elif file.endswith('json'):
            with open(file, 'r') as json_file:
                to_return = json.load(json_file)

        elif file.endswith('csv'):
            with open(file, 'r') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for idx, section in enumerate(csv_reader):
                    if section[4] in ['L-bulb', 'bulb', 'hp']:
                        to_return[str(idx)] = {'stf_web_height': [float(section[0]) - float(section[3]), 'm'],
                                               'stf_web_thk': [float(section[1]),'m'],
                                               'stf_flange_width': [float(section[2]),'m'],
                                               'stf_flange_thk': [float(section[3]), 'm'],
                                               'stf_type': [section[4], '']}
                    else:
                        to_return[str(idx)] = {'stf_web_height': [float(section[0]), 'm'],
                                               'stf_web_thk': [float(section[1]),'m'],
                                               'stf_flange_width': [float(section[2]),'m'],
                                               'stf_flange_thk': [float(section[3]), 'm'],
                                               'stf_type': [section[4], '']}

    if to_json:
        with open('sections.json', 'w') as file:
            json.dump(to_return, file)

    if to_csv is not None:
        with open(to_csv, 'w', newline = '') as file:
            section_writer = csv.writer(file)
            for line in return_csv:
                section_writer.writerow(line)
    if obj is not None:  # This will return a modified object.
        if type(obj) is not list:
            obj = [obj, ]
            append_list = [[],]
        else:
            append_list = [list() for dummy in obj]
    else:
        append_list = list()

    for key, value in to_return.items():
        if obj is not None:  # This will return a modified object.
            for idx, iter_obj in enumerate(obj):
                new_obj = copy.deepcopy(iter_obj)
                new_obj_prop = new_obj.get_structure_prop()
                for prop_name, prop_val in value.items():
                    new_obj_prop[prop_name] = prop_val
                new_obj.set_main_properties(new_obj_prop)
                append_list[idx].append(new_obj)
        else:
            to_return_final.append(value)
    if len(append_list) == 1:
        to_return_final = append_list[0]
    elif len(append_list) == 0:
        pass
    elif len(append_list) > 1:
        to_return_final = append_list

    return to_return_final
