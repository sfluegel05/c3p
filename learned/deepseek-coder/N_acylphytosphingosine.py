"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the phytosphingosine backbone pattern
    # Phytosphingosine typically has a long carbon chain with multiple hydroxyl groups and an amide bond
    phytosphingosine_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OH])[CHX4]([OH])[CH2X4][NX3][CX3](=[OX1])")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Look for the fatty acyl group attached to the nitrogen via an amide bond
    fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CH2X4]")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl group attached to nitrogen found"

    # Check for the presence of hydroxyl groups on the sphingosine backbone
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 3"

    # Check molecular weight - N-acylphytosphingosine typically has a high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for N-acylphytosphingosine"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for N-acylphytosphingosine"
    if o_count < 4:
        return False, "Must have at least 4 oxygens (hydroxyl groups and amide bond)"

    return True, "Contains phytosphingosine backbone with a fatty acyl group attached to the nitrogen and hydroxyl groups"