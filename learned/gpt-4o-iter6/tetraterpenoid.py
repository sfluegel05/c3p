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
    
    # Count carbon atoms - setting a more specific range for tetraterpenoids (38 to 42 is typical)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 38 or carbon_count > 42:
        return False, "Carbon count out of tetraterpenoid range"
    
    # Look for long conjugated carbon chains with specific types of bonds or even rings
    # This is a simplistic pattern for detecting polyconjugation, focusing on C=C double bonds
    long_chain_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long conjugated chain pattern typical for tetraterpenoids found"

    # Besides simple ketone or hydroxyl, check for presence of epoxides and glycosides 
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    hydroxyl_pattern = Chem.MolFromSmarts("O[CH]")
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]([C@H])")
    
    if (mol.HasSubstructMatch(ketone_pattern) or 
        mol.HasSubstructMatch(hydroxyl_pattern) or 
        mol.HasSubstructMatch(epoxide_pattern) or
        mol.HasSubstructMatch(glycoside_pattern)):
        return True, "Contains tetraterpenoid-like C40 structure with notable functional modifications"

    # Fallback if ambiguous
    return False, "Structure does not conclusively match typical tetraterpenoid characteristics"

# Example usage
is_tetraterpenoid("CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C")  # gamma-carotene example