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

    # Phosphate group detection: P(=O)(O)(O) with potential external ester linkage
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group incorrectly configured"

    # Ester group linked to glycerol: [C]-(C=O)-O with [C] in the scaffold of glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]1COC1")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages with glycerol, need exactly 2"

    # Check glycerol backbone (checked in context of esters)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO)OC1C=O") # Simplified part of glycerol to ensure checking esterification context
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing glycerol backbone in expected configuration"

    # Verify correct attachment of phosphoric acid to glycerol (ensure linkage)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15: # Phosphorus
            num_oxygen_bonds = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8) # Count oxygen bonds
            if num_oxygen_bonds != 3: # Expect exactly 3 oxygen atoms attached to the phosphate
                return False, "Phosphate group has incorrect oxygen attachments"

    return True, "Structure matches phosphatidic acid with glycerol backbone and appropriate esterifications"