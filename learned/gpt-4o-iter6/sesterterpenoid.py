"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are derived from sesterterpenes and can have a modified C25 skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms, consider modified C25 skeletons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 45:
        return False, f"Carbon count ({c_count}) not typical or too large for a sesterterpenoid"

    # Use SMARTS to define key fragments or motifs - multi-ring systems typical in sesterterpenoids
    complex_cyclic_patterns = [
        Chem.MolFromSmarts("C1CCC2CC3CCC(C1)C23"),  # Three-fused rings
        Chem.MolFromSmarts("C1C=C2CCC3(C)C=C(C)CCC3(C)CC2C1"),  # Sesterterpenoid-like
    ]
    
    # Check for presence of complex cyclic structure typical in sesterterpenoids
    if not any(mol.HasSubstructMatch(pattern) for pattern in complex_cyclic_patterns):
        return False, "No key cyclic structures typical of sesterterpenoids found"
    
    # Calculate and validate additional characteristics like molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1000:  # Sesterterpenoids generally within this range
        return False, "Molecular weight not within typical range for sesterterpenoids"

    # If structure fits broad criteria
    return True, "Likely a sesterterpenoid based on carbon count and presence of key cyclic structures"