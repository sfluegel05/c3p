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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 35:  # Adjusted range to account for variations
        return False, f"Carbon count ({c_count}) not within typical range for sesterterpenoid"

    # Expanded set of complex cyclic patterns typical in sesterterpenoids
    complex_cyclic_patterns = [
        Chem.MolFromSmarts("C1CCC2CCC(C1)C2"),  # Simple bicyclic system
        Chem.MolFromSmarts("C1CCCC2CCCC3(C1)C2CCC3"),  # Tricyclic system
        Chem.MolFromSmarts("C1CCC2CC3CC(C2C1)CC3"),  # Larger tricyclic system
        Chem.MolFromSmarts("C1C[C@H]2CC[C@@H]3C[C@H](C1)C2(C)CC3"),  # ABC-ring terpenoid scaffold
    ]

    # Check for presence of any complex cyclic structures typical in sesterterpenoids
    if not any(mol.HasSubstructMatch(pattern) for pattern in complex_cyclic_patterns):
        return False, "No key cyclic structures typical of sesterterpenoids found"
    
    # Calculate and validate additional characteristics like molecular weight potentially
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 700:  # Adjusted weight range based on typical sesterterpenoids
        return False, "Molecular weight not within typical range for sesterterpenoids"

    # Extra check for typical functional groups, if possible, can be added.

    # If structure fits broad criteria
    return True, "Contains key characteristics of sesterterpenoids"