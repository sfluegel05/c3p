"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check number of carbon atoms, usually tetraterpenoids have around 40
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 35:
        return False, f"Too few carbon atoms for tetraterpenoid: {num_carbons} carbons"

    # Check for the presence of unsaturated bonds, characteristic of terpenes
    num_double_bonds = Descriptors.CalcNumAromaticBonds(mol) + Descriptors.CalcNumAliphaticDoubleBonds(mol)
    if num_double_bonds < 10:  # arbitrary cutoff for a complex terpene structure
        return False, "Insufficient double bonds for a tetraterpenoid structure"
    
    # Scan for specific terpenoid functional groups or motifs (indicative, not exhaustive)
    # Example: Alcohol groups, indicative of oxygenation
    alcohol_pattern = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Lacks common terpenoid modification groups (e.g., alcohol)"
    
    # Pattern matching for characteristic skeleton modification or restructuring could be complex
    # Here we simply return None for structures that require in-depth analysis beyond simple rule matching

    # If initial checks pass, assume tetraterpenoid unless specific exclusions apply
    return True, "Structure likely contains the characteristic features of a tetraterpenoid"

# Example use: is_tetraterpenoid("CC(=O)NC1=CC=C(O)C=C1")