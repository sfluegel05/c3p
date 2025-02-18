"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is defined as any organic compound with an initial boiling point ≤250°C at standard atmospheric pressure.
    Uses molecular weight and functional groups as heuristic indicators.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check if organic (contains carbon)
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if not has_carbon:
        return False, "Not an organic compound"
    
    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)
    
    # Define functional groups that increase boiling point
    high_bp_groups = [
        ('hydroxyl', '[OH]'),
        ('carboxylic_acid', '[CX3](=O)[OX2H1]'),
        ('amine_primary', '[NH2]'),
        ('nitro', '[N+](=O)[O-]'),
        ('ester', '[OX2][CX3](=[OX1])[#6]'),
    ]
    
    # Check for presence of high boiling point groups
    has_high_bp = False
    for name, smarts in high_bp_groups:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_high_bp = True
            break
    
    # Apply heuristic rules
    if mol_wt > 300:
        return False, f"Molecular weight {mol_wt:.2f} > 300"
    elif has_high_bp and mol_wt > 200:
        return False, f"Has {name} group and molecular weight {mol_wt:.2f} > 200"
    else:
        return True, f"Molecular weight {mol_wt:.2f} and functional groups indicate boiling point ≤250°C"