"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are characterized by a flavan-3-ol backbone and substituted hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined SMARTS pattern for flavan-3-ol backbone: C6-C3-C6 pattern with hydroxyl on ring C 
    flavan_3_ol_pattern = Chem.MolFromSmarts("c1cc2OC(C(O)CC2c1)c3ccc(O)cc3")
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol backbone found"
    
    # Check additional substituents
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    methoxy_pattern = Chem.MolFromSmarts("CO")
    ester_pattern = Chem.MolFromSmarts("O=C(O)")

    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    num_oxygen_containing_subs = len(hydroxyl_matches) + len(methoxy_matches) + len(ester_matches)
    
    if num_oxygen_containing_subs < 3:
        return False, f"Insufficient oxygen-containing groups, found {num_oxygen_containing_subs}"

    # Check stereochemistry: Catechins should have certain chiral centers
    expected_chiral_centers = [(2, 'R'), (3, 'S')]
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    valid_chirality = all((center in chiral_centers) for center in expected_chiral_centers)
    if not valid_chirality:
        return False, "Chiral centers do not match expected catechin configuration"

    return True, "Contains flavan-3-ol backbone with sufficient oxygen-containing groups and correct stereochemistry"