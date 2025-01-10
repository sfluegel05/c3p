"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Count carbon atoms to match known examples of tetraterpenoid structures
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (30 <= carbon_count <= 45):
        return False, "Carbon count out of typical tetraterpenoid range"

    # Check for long conjugated polyene or similar patterns
    polyene_patterns = [
        Chem.MolFromSmarts("C=C-C=C"),  # Simple diene unit
    ]
    if not any(mol.HasSubstructMatch(patt) for patt in polyene_patterns):
        return False, "No elongated polyene chain pattern typical for tetraterpenoids found"
        
    # Additional functional modification checks
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    epoxide_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    glycoside_patterns = [
        Chem.MolFromSmarts("O[C@H]"),  # Basic sugar attachment
        Chem.MolFromSmarts("[C@H]1O[C@H]"),  # More comprehensive sugar variants
    ]

    # Need at least one of the below common functionality to tag as tetraterpenoid
    if any(mol.HasSubstructMatch(patt) for patt in [ketone_pattern, hydroxyl_pattern, epoxide_pattern] + glycoside_patterns):
        return True, "Contains tetraterpenoid-like C40 structure with notable functional modifications"

    return False, "Structure does not conclusively match typical tetraterpenoid characteristics"

# Example usage
print(is_tetraterpenoid("CC(C)CCC\\C(C)=C\\CC\\C(C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(/C)\\C=C\\C=C(/C)\\C=C\\C=C(/C)CCCC(C)(C)O"))  # Chloroxanthin for review