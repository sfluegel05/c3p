"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde with a long aliphatic carbon chain derived from fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aldehyde group [CX3H]=O at the end
    # Aldehyde group pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No terminal aldehyde group found"

    # Check for a long carbon chain
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 6:
        return False, f"Carbon chain too short for fatty aldehyde (found {num_carbons} carbons)"

    # Check molecular weight to ensure it's consistent with fatty aldehydes
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:  # Rough estimate, could be adjusted based on chain length
        return False, "Molecular weight too low for a typical fatty aldehyde"
    
    # Ensure there are no conflicting functional groups that invalidate the molecule.
    # WARNING: Simplified check. More comprehensive checks depend on other aldehyde inhibitors.
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][#8]")):
        return False, "Unexpected functional groups present"

    return True, "Valid fatty aldehyde: Terminal aldehyde group with a long aliphatic chain"