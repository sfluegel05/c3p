"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone consisting of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for butenolide substructure
    butenolide_pattern = Chem.MolFromSmarts("[O-,O]1C=CC(=O)O1")
    if not mol.HasSubstructMatch(butenolide_pattern):
        return False, "No butenolide substructure found"
    
    # Check molecular weight range (typically 80-300 Da for butenolides)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 300:
        return False, f"Molecular weight ({mol_wt:.2f} Da) out of typical range for butenolides"
    
    # Check elemental composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 4 or o_count < 2:
        return False, "Insufficient carbon or oxygen atoms for butenolide"
    
    return True, "Contains a 2-furanone skeleton with a gamma-lactone arrangement (butenolide)"