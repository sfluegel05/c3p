"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Minerals as defined by geological processes, including ionic compounds, covalent structures, and elemental forms.
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Minerals are typically:
    - Ionic compounds with metal cations and inorganic anions
    - Covalent structures of metals with nonmetals (O, S, F, Cl, etc.)
    - Elemental forms of metals/metalloids or certain nonmetals (e.g., sulfur)
    - Hydrates where the core compound meets mineral criteria

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define metals/metalloids and allowed nonmetals for covalent structures
    metals = {'Li','Be','Na','Mg','Al','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',
              'Ni','Cu','Zn','Ga','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
              'Ag','Cd','In','Sn','Sb','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
              'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir',
              'Pt','Au','Hg','Tl','Pb','Bi','Po','Fr','Ra','Ac','Th','Pa','U','Si','B','As','Se','Te','Ge'}
    covalent_nonmetals = {'O','S','F','Cl','Br','I','N','P'}

    # Split into fragments (separate ions/metal complexes)
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    
    # Remove water fragments for hydrate handling
    non_water_frags = []
    for frag in fragments:
        formula = Chem.rdMolDescriptors.CalcMolFormula(frag)
        if formula == 'H2O':
            continue
        non_water_frags.append(frag)
    
    # Check 1: Ionic compounds (multiple fragments with metal cations and inorganic anions)
    if len(non_water_frags) >= 2:
        has_metal_cation = False
        has_inorganic_anion = False
        
        for frag in non_water_frags:
            frag_charge = sum(atom.GetFormalCharge() for atom in frag.GetAtoms())
            contains_metal = any(atom.GetSymbol() in metals for atom in frag.GetAtoms())
            contains_carbon = any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms())
            
            if frag_charge > 0 and contains_metal:
                has_metal_cation = True
            elif frag_charge < 0 and not contains_carbon:  # Anion must be inorganic (no carbon)
                has_inorganic_anion = True
        
        if has_metal_cation and has_inorganic_anion:
            return True, "Ionic compound with metal cation and inorganic anion"

    # Check 2: Covalent structures (single fragment with metal + nonmetal bonds)
    if len(non_water_frags) == 1:
        has_metal = False
        has_nonmetal = False
        for atom in non_water_frags[0].GetAtoms():
            sym = atom.GetSymbol()
            if sym in metals:
                has_metal = True
            elif sym in covalent_nonmetals:
                has_nonmetal = True
        
        if has_metal and has_nonmetal:
            return True, "Covalent metal-nonmetal compound"

    # Check 3: Elemental forms (single element, neutral)
    elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if len(elements) == 1:
        elem = elements.pop()
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        if (elem in metals or elem in {'S','C'}) and total_charge == 0:
            return True, "Elemental form of a mineral"

    return False, "Does not meet mineral criteria"