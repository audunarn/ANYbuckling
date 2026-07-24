'''
Runner for the licensed DNV PULS (Panel Ultimate Limit State) Excel
sheet. Requires the ``excel`` extra (xlwings) and a PULS Excel sheet to
actually execute runs; importing this module does not require xlwings.
'''
from .panel import PULSpanel

__all__ = ["PULSpanel"]
