"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes, which have a modified C40 carbon skeleton often with functional modifications.
    
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
    
    # Count carbon atoms in the range expected for tetraterpenoid base structures
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Checking for base carbon backbone - note leniency for modifications
    if not (36 <= carbon_count <= 44):
        return False, "Carbon count out of typical tetraterpenoid range"

    # Improved pattern checking for long conjugated systems
    # Multiple conjugated trans double bonds (general form)
    polyene_chain_pattern = Chem.MolFromSmarts("C=C-C=C-C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyene_chain_pattern):
        return False, "No elongated polyene chain pattern typical for tetraterpenoids found"

    # Check for characteristic functional modifications
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    epoxide_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]")  # More general pattern for sugars

    # Need at least one of the below common functionality to tag as tetraterpenoid
    if any(mol.HasSubstructMatch(patt) for patt in [ketone_pattern, hydroxyl_pattern, epoxide_pattern, glycoside_pattern]):
        return True, "Contains tetraterpenoid-like C40 structure with notable functional modifications"

    return False, "Structure does not conclusively match typical tetraterpenoid characteristics"

# Example usage
print(is_tetraterpenoid("CC(C)CCC\\C(C)=C\\CC\\C(C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(/C)\\C=C\\C=C(/C)\\C=C\\C=C(/C)CCCC(C)(C)O"))  # Chloroxanthin for review