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

    # Look for phytosphingosine backbone pattern
    # [CH2OH]-[CH(NH)]-[CH(OH)]-[CH(OH)]-[CH2]- connected to long chain
    phytosphingosine_pattern = Chem.MolFromSmarts("[CH2][OH].[CH]([NH])[CH]([OH])[CH]([OH])[CH2]")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Count hydroxyl groups (should have at least 3)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 3:
        return False, f"Insufficient hydroxyl groups (found {hydroxyl_matches}, need ≥3)"

    # Check for long carbon chain attached to amide
    # First, get the amide carbons
    amide_carbons = mol.GetSubstructMatches(amide_pattern)
    if not amide_carbons:
        return False, "No amide group carbon found"
    
    # Count carbons in molecule (should have enough for backbone and acyl chain)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Minimum carbons for basic structure
        return False, f"Insufficient carbon atoms for N-acylphytosphingosine (found {c_count}, need ≥20)"

    # Check for aliphatic chain characteristic of fatty acyl group
    fatty_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_chain_pattern):
        return False, "No fatty acyl chain found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Insufficient rotatable bonds for N-acylphytosphingosine structure"

    return True, "Contains phytosphingosine backbone with N-acyl group and required hydroxyl groups"