"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes, which typically have a modified C40 carbon skeleton and various functionalities.
    
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
    
    # Count carbon atoms - adjusted range for potential rearrangements/modifications
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 36 or carbon_count > 44:
        return False, "Carbon count out of typical tetraterpenoid range"
    
    # Look for terpenoid-like conjugated polyene chains (revised pattern)
    polyene_chain_pattern = Chem.MolFromSmarts("[CH]=[CH]-[CH]=[CH]-[CH]=[CH]-[CH]=[CH]")
    if not mol.HasSubstructMatch(polyene_chain_pattern):
        return False, "No elongated polyene chain pattern typical for tetraterpenoids found"

    # Consider the presence of functional groups: ketones, hydroxyls, epoxides, glycosides
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    epoxide_pattern = Chem.MolFromSmarts("[OX2]1[CX4][OX2]1")
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]([C@H])")

    # Match at least one of the functional pattern or a combination
    if (mol.HasSubstructMatch(ketone_pattern) or
        mol.HasSubstructMatch(hydroxyl_pattern) or
        mol.HasSubstructMatch(epoxide_pattern) or
        mol.HasSubstructMatch(glycoside_pattern)):
        return True, "Contains tetraterpenoid-like C40 structure with notable functional modifications"

    # Fallback if ambiguous
    return False, "Structure does not conclusively match typical tetraterpenoid characteristics"

# Example usage
is_tetraterpenoid("CC(C)CCC\\C(C)=C\\CC\\C(C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(/C)\\C=C\\C=C(/C)\\C=C\\C=C(/C)CCCC(C)(C)O")  # Review of example chloroxanthin