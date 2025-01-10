"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes, which have a C40 carbon skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms - tetraterpenoids have a large number of carbons, theoretically around 40
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_count < 35:
        return False, "Too few carbons for a tetraterpenoid-like structure"

    # Look for long conjugated carbon chain - typical in tetraterpenoids
    # A very simple pattern for long carbon chains could be represented with repeating [CH]
    long_chain_pattern = Chem.MolFromSmarts("[*]~[*]~[*]~[*]~[*]~[*]~[*]~[*]~[*]~[*]") # Simplistic long chain
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long conjugated carbon chain found"

    # Additionally, look for some specific features like ketones or hydroxyls
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    hydroxyl_pattern = Chem.MolFromSmarts("O[CH]")
    if mol.HasSubstructMatch(ketone_pattern) or mol.HasSubstructMatch(hydroxyl_pattern):
        return True, "Contains tetraterpenoid-like C40 structure with functional modifications"

    # Fallback if ambiguous
    return False, "Structure does not match typical tetraterpenoid characteristics"

# Example usage
is_tetraterpenoid("CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C")  # gamma-carotene example