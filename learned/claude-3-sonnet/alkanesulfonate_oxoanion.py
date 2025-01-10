"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a sulfonate group (-SO3-) attached to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonate group pattern: -S([O-])(=O)=O
    sulfonate_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])[O-]")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "No sulfonate group found"

    # Get sulfonate group matches
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    
    # Check each sulfonate group
    valid_sulfonate = False
    for match in sulfonate_matches:
        sulfur_idx = match[0]  # Get sulfur atom index
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Check neighbors of sulfur
        for neighbor in sulfur_atom.GetNeighbors():
            # If neighbor is carbon, we have an alkanesulfonate
            if neighbor.GetAtomicNum() == 6:  # Carbon atomic number
                valid_sulfonate = True
                break
                
    if not valid_sulfonate:
        return False, "Sulfonate group not attached to carbon"

    # Count formal charges to ensure single negative charge on sulfonate
    total_negative_charge = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    if total_negative_charge < 1:
        return False, "No negative charge found"

    # Additional validation: check for S-C bond explicitly
    sc_bond_pattern = Chem.MolFromSmarts("[#6]-[SX4](=[OX1])(=[OX1])[O-]")
    if not mol.HasSubstructMatch(sc_bond_pattern):
        return False, "No carbon-sulfur bond found in sulfonate group"

    return True, "Contains sulfonate group (-SO3-) attached to carbon"