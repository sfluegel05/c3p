"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:26873 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import Counter

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is a terpenoid derived from a sesquiterpene (C15 hydrocarbons built from three isoprene units).
    The skeleton may be rearranged or modified by the removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get elemental composition
    atoms = mol.GetAtoms()
    elem_counts = Counter()
    for atom in atoms:
        elem_counts[atom.GetSymbol()] += 1
    C = elem_counts.get('C', 0)
    H = elem_counts.get('H', 0)
    N = elem_counts.get('N', 0)
    O = elem_counts.get('O', 0)
    halogens = elem_counts.get('F',0)+elem_counts.get('Cl',0)+elem_counts.get('Br',0)+elem_counts.get('I',0)
    
    # Sesquiterpenoids typically have around 15 carbons, allow some variation
    if C < 13 or C > 17:
        return False, f"Carbon count is {C}, which is not typical for sesquiterpenoids"
    
    # Calculate the degree of unsaturation
    unsaturation = (2 * C + 2 + N - H - halogens) / 2
    if unsaturation < 4:
        return False, f"Degree of unsaturation is {unsaturation}, which is low for sesquiterpenoids"

    # Check for isoprene units or their rearranged forms
    # Define SMARTS patterns for isoprene units and common sesquiterpene skeletons
    isoprene_smarts_list = [
        'C(=C)C-C=C',           # Standard isoprene unit
        'C-C(=C)-C=C',          # Rearranged isoprene
        'C=C-C=C-C',            # Another rearranged isoprene
        'C=C-C-C-C=C',          # Extended isoprene pattern
    ]
    found_isoprene = False
    for smarts in isoprene_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            found_isoprene = True
            break
    if not found_isoprene:
        return False, "No isoprene units or sesquiterpene skeletons detected"

    # Check for common sesquiterpene skeletons
    sesquiterpene_skeletons = [
        'C1CCC2(C1)CCC(C=C2)C',  # Germacrane skeleton
        'C1=CC=CC2(C1)CCCCC2',   # Cadinane skeleton
        'C1CCC2=C(C1)CCCC2',     # Eudesmane skeleton
        'C1CC2CC3CC(C3)CC2C1',   # Guaiane skeleton
        'C1CC2CCC3=C(CC3)C2C1',  # Aromadendrane skeleton
    ]
    found_skeleton = False
    for smarts in sesquiterpene_skeletons:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            found_skeleton = True
            break
    if not found_skeleton:
        return False, "Molecule does not match common sesquiterpene skeletons"

    # If the molecule passes all checks
    return True, "Molecule meets criteria for a sesquiterpenoid"