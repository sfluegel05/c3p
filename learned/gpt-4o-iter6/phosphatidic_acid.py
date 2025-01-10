"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is characterized by a glycerol backbone where one hydroxy group,
    commonly but not necessarily primary, is esterified with phosphoric acid,
    and the other two are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone with hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("OC[C@@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing glycerol backbone with hydroxyls"

    # Phosphate group linking pattern: P(=O)(O)O
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, "Phosphate group incorrectly configured"

    # Check for two ester linkages: C(=O)O attached to glycerol
    ester_pattern = Chem.MolFromSmarts("[C][O][C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 2 for fatty acids"

    # Verify no additional linkages to the phosphate that would indicate another group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15: # Phosphorus
            num_oxygen_bonds = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8) # Count oxygen bonds
            if num_oxygen_bonds != 3: # Expect exactly 3 oxygen atoms attached to the phosphate
                return False, "Phosphate group has incorrect oxygen attachments"

    return True, "Structure matches phosphatidic acid with glycerol backbone and appropriate esterifications"