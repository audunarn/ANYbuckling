'''
Packaging-level tests: imports, exports, packaged data and CLI wiring.
'''
import os


def test_top_level_exports():
    import anybuckling

    for name in anybuckling.__all__:
        assert getattr(anybuckling, name, None) is not None, name
    assert anybuckling.__version__


def test_bulb_section_data_is_packaged():
    from anybuckling.sections import bulb_section_file, helper_read_section_file

    path = bulb_section_file()
    assert os.path.isfile(path)
    sections = helper_read_section_file(path)
    assert len(sections) > 100
    assert {'stf_web_height', 'stf_web_thk', 'stf_flange_width',
            'stf_flange_thk', 'stf_type'} <= set(sections[0])


def test_puls_module_imports_without_xlwings():
    from anybuckling.puls import PULSpanel

    runner = PULSpanel()
    assert runner.puls_acceptance == 0.87


def test_semianalytical_cli_parser():
    from anybuckling.semianalytical.solver import _build_argument_parser

    parser = _build_argument_parser()
    args = parser.parse_args(['benchmark', '--csv', 'dummy.csv'])
    assert args.csv == 'dummy.csv'
