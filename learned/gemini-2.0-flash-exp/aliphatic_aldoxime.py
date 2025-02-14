"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is an aldoxime where the carbon attached to the oxime group is derived from an aliphatic aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check for the oxime group. The carbon needs to have a single hydrogen.
    oxime_pattern = Chem.MolFromSmarts("[CH1]=[N]-O")
    if not mol.HasSubstructMatch(oxime_pattern):
        return False, "No oxime group found"
        
    # Get the carbon atoms from the oxime group
    matches = mol.GetSubstructMatches(oxime_pattern)
    
    # 3. Check for aliphatic nature by checking that no neighboring atoms are part of aromatic ring
    for match in matches:
        carbon_index = match[0]  # the carbon index is always the first in the match tuple
        carbon_atom = mol.GetAtomWithIdx(carbon_index)
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
              return False, "The carbon of the oxime group is attached to aromatic ring"

    # 4. Return True if all conditions are met
    return True, "Molecule is an aliphatic aldoxime"