"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime contains an aliphatic chain and an aldehyde group
    converted to an oxime (R-C=NOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for the classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the oxime group pattern "-C=NO"
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OH]")
    if not mol.HasSubstructMatch(oxime_pattern):
        return False, "No oxime group found"

    # Ensure the carbon in C=NO is part of a non-aromatic chain, i.e., aliphatic
    for match in mol.GetSubstructMatches(oxime_pattern):
        carbon_atom = match[0]
        
        # Check if the carbon is part of an aliphatic chain
        atom = mol.GetAtomWithIdx(carbon_atom)
        if atom.GetIsAromatic():
            return False, "Carbon in C=NO is part of an aromatic ring"

    return True, "Contains an aliphatic oxime group"