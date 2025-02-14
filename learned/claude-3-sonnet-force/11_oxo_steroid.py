"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: CHEBI:35665 11-oxo steroid
Any oxo steroid that has an oxo substituent at position 11.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone pattern
    # Use a more general pattern to cover diverse steroid structures
    steroid_pattern = Chem.MolFromSmarts("[C@]12CCC[C@]3([C@@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for oxo group at position 11
    # Use a more flexible pattern to capture different ring conformations
    oxo_pattern = Chem.MolFromSmarts("[C@]12CCC[C@@]3([C@@H]([C@@H]1CC[C@@H]2C(=O))CCC4=CC(=O)CC[C@]34C)C")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group found at position 11"

    # Check molecular weight range (250 - 600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.2f} is out of range for 11-oxo steroids"

    # Count carbons (20-30) and oxygens (2-5)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20 or c_count > 30:
        return False, f"Carbon count {c_count} is out of range for 11-oxo steroids"
    if o_count < 2 or o_count > 5:
        return False, f"Oxygen count {o_count} is out of range for 11-oxo steroids"

    return True, "Contains a steroid backbone with an oxo group at position 11"