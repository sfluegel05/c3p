"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:28789 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are organic compounds whose structure is based on derivatives of a phenyl-substituted 1-phenylpropane
    possessing a C6-C3-C6 skeleton or a related structure, often with additional rings, substituents, or oligomers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavonoid backbone patterns
    flavonoid_patterns = [
        Chem.MolFromSmarts("[c1c-,c;c,c]1[C@@H](C)[C@H](C)[c2c-,c;c,c]"), # C6-C3-C6 core
        Chem.MolFromSmarts("[c1c-,c;c,c]1[C@H](C)[C@@H](C)[c2c-,c;c,c]"), # C6-C3-C6 isomer
        Chem.MolFromSmarts("[c1c-,c;c,c]1[C@@H](C)[C@@H](C)[c2cc-,ccc2]"), # C6-C3-C3 core
        Chem.MolFromSmarts("[c1c-,c;c,c]1[C@@H](C)[C@H](C)[c2cc-,ccc2]"), # C6-C3-C3 isomer
        Chem.MolFromSmarts("[c1c-,c;c,c]1[C@@H](C)[C@@H](C)[c2c-,c;c,c][c3c-,c;c,c]"), # C6-C3-C6-C3
        Chem.MolFromSmarts("[c1c-,c;c,c]1[C@@H](C)[C@H](C)[c2c-,c;c,c][c3c-,c;c,c]") # C6-C3-C6-C3 isomer
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_patterns):
        return False, "No flavonoid backbone found"
    
    # Look for common flavonoid substituents
    substituents = [
        Chem.MolFromSmarts("[OX2H]"),  # Hydroxyl
        Chem.MolFromSmarts("[OX2CH3]"),  # Methoxy
        Chem.MolFromSmarts("[C=C](C)C"),  # Prenyl
        Chem.MolFromSmarts("[OX2C1OCC(O)C(O)C1]"),  # Sugar moiety
    ]
    
    if not any(mol.HasSubstructMatch(sub) for sub in substituents):
        return False, "No common flavonoid substituents found"
    
    # Check for aromatic rings
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if ring.isAromatic()]
    if len(aromatic_rings) < 2:
        return False, "Not enough aromatic rings for flavonoid"
    
    # Check molecular weight - flavonoids typically 200-800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.2f} outside typical range for flavonoids"
    
    return True, "Contains flavonoid backbone and common substituents"