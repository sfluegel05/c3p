"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Minerals as defined by geological processes, including ionic compounds and certain covalent structures.
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Minerals are typically ionic compounds with metal cations and inorganic anions,
    or covalent structures like oxides/sulfides of metals/metalloids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of elements considered metals/metalloids for mineral classification
    metals = {'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
              'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
              'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
              'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
              'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Si', 'B', 'As', 'Se', 'Te', 'Ge'}
    
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    has_cation = False
    has_anion = False
    
    # Check for ionic compounds (multiple fragments with charges)
    for frag in fragments:
        for atom in frag.GetAtoms():
            symbol = atom.GetSymbol()
            charge = atom.GetFormalCharge()
            # Check for metal cations
            if symbol in metals and charge > 0:
                has_cation = True
            # Check for anions
            if charge < 0:
                has_anion = True
                
    if has_cation and has_anion:
        return True, "Ionic compound with metal cation and anion"
    
    # Check for covalent metal oxides/sulfides (single fragment)
    if len(fragments) == 1:
        metal_present = False
        oxygen_count = 0
        sulfur_count = 0
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in metals:
                metal_present = True
            if symbol == 'O':
                oxygen_count += 1
            if symbol == 'S':
                sulfur_count += 1
        
        # Metal oxide (e.g., SiO2, Fe2O3)
        if metal_present and oxygen_count > 0:
            return True, "Covalent metal oxide"
        # Metal sulfide (e.g., FeS2 as covalent structure)
        if metal_present and sulfur_count > 0:
            return True, "Covalent metal sulfide"
    
    # Check for elemental forms (e.g., native metals, sulfur)
    elements_present = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if len(elements_present) == 1:
        elem = elements_present.pop()
        if elem in metals or elem in {'S', 'C'}:
            return True, "Elemental form of a mineral"
    
    return False, "Does not meet mineral criteria"