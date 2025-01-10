"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_N_acylphytosphingosine, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatches(amide_pattern):
        return False, "No amide group found"

    # Look for phytosphingosine core with specific stereochemistry
    # 2S,3S,4R stereochemistry is required for phytosphingosine
    phyto_core = Chem.MolFromSmarts("[CH2][OX2H1]-[CH]([NX3])-[CH]([OX2H1])-[CH]([OX2H1])-[CH2]")
    if not mol.HasSubstructMatch(phyto_core):
        return False, "No phytosphingosine core structure found"

    # Count primary and secondary hydroxyl groups
    primary_oh_pattern = Chem.MolFromSmarts("[CH2][OX2H1]")
    secondary_oh_pattern = Chem.MolFromSmarts("[CH]([#6])[OX2H1]")
    
    primary_oh_count = len(mol.GetSubstructMatches(primary_oh_pattern))
    secondary_oh_count = len(mol.GetSubstructMatches(secondary_oh_pattern))
    
    if primary_oh_count < 1 or secondary_oh_count < 2:
        return False, f"Incorrect hydroxyl pattern (need 1 primary, 2 secondary OH groups)"

    # Check for long carbon chain characteristic of fatty acyl group
    fatty_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No fatty acyl chain found"

    # Count carbons in molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:  # Minimum carbons for basic structure
        return False, f"Insufficient carbon atoms (found {c_count}, need ≥18)"

    # Check for characteristic N-acyl linkage pattern
    n_acyl_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]-[CX4]")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl linkage pattern found"

    # Allow for sugar modifications by checking if extra oxygens are present in cyclic form
    sugar_pattern = Chem.MolFromSmarts("[C]1[O][C]([C])[C][C]([O])[C]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    
    # Count rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8 and not has_sugar:
        return False, "Insufficient chain length"

    return True, "Contains phytosphingosine backbone with N-acyl group and required hydroxyl groups"