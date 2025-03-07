import numpy as np
import requests
import sys

# ----------------------------------------------------------------------------------------------------
def cross_section(A, Z):

    elements = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 
        11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 
        20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 
        29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 
        38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 
        47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 
        56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 
        65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 
        74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb'
    }    
    
    element = elements.get(Z)

    if not element:
        print(f'Element with atomic number {Z} not found.')
        return np.array([]), np.array([])

    url = f'https://tendl.web.psi.ch/tendl_2023/gamma_file/{element}/{element}{A:03}/tables/xs/nonelastic.tot'

    response = requests.get(url)

    if response.status_code != 200:
        print(f'Failed to retrieve data: {response.status_code}')
        return np.array([]), np.array([])

    data_lines = response.text.strip().split('\n')

    eps = []
    xs = []

    for line in data_lines:
        if line.startswith('#'):
            continue
        
        try:
            values = list(map(float, line.split()))
            if len(values) >= 2:
                eps.append(values[0])
                xs.append(values[1])
        except ValueError:
            print(f'Invalid data found in line: {line}')

    return np.array(eps), np.array(xs)

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    
    A = int(sys.argv[1])
    Z = int(sys.argv[2])

    eps, xs = cross_section(A, Z)

    print(eps.tolist())  
    print(xs.tolist())

# ----------------------------------------------------------------------------------------------------
