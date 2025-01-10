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
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Look for phytosphingosine core 
    # [CH2OH]-[CH](NH)-[CH](OH)-[CH](OH)-[CH2]-
    phyto_core = Chem.MolFromSmarts("[CH2][OX2]-[CH]([NX3])-[CH]([OX2])-[CH]([OX2])-[CH2]")
    if not mol.HasSubstructMatch(phyto_core):
        return False, "No phytosphingosine core structure found"

    # Count hydroxyl groups
    primary_oh_pattern = Chem.MolFromSmarts("[CH2][OX2H1]")
    secondary_oh_pattern = Chem.MolFromSmarts("[CH]([#6])[OX2H1]")
    
    if not mol.HasSubstructMatch(primary_oh_pattern):
        return False, "Missing primary hydroxyl group"
    
    secondary_matches = len(mol.GetSubstructMatches(secondary_oh_pattern))
    if secondary_matches < 2:
        return False, f"Insufficient secondary hydroxyl groups (found {secondary_matches}, need ≥2)"

    # Check for fatty acyl chain
    # Look for at least 6 carbons in a chain attached to the amide
    fatty_chain = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No suitable fatty acyl chain found"

    # Count carbons and check molecular mass
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Insufficient carbon atoms (found {c_count}, need ≥18)"

    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 300:  # Minimum weight for basic structure
        return False, f"Molecular weight too low ({mol_weight:.1f} Da)"

    # Check for possible sugar modifications
    sugar_pattern = Chem.MolFromSmarts("[C]1[O][C]([C])[C][C]([O])[C]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    
    # Additional checks for chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8 and not has_sugar:
        return False, "Chain length too short for N-acylphytosphingosine"

    # Check nitrogen count (should typically have 1, but allow for modifications)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen found"
    if n_count > 3 and not has_sugar:  # Allow more N if sugar modifications present
        return False, f"Too many nitrogens for basic structure (found {n_count})"

    return True, "Contains phytosphingosine backbone with N-acyl group and required hydroxyl groups"