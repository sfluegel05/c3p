"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    Proanthocyanidins are characterized by flavan-3-ol units, aromatic hydroxyl groups, 
    and specific linkages between these units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define structural search patterns
    # Flavan-3-ol core structure
    flavan3ol_patterns = [
        Chem.MolFromSmarts('Oc1cc(O)c2C[C@@H]([C@H](O)C2Oc2cc(O)cc(O)c12)C'),
        Chem.MolFromSmarts('Oc1cc(O)c2C[C@H]([C@@H](O)C2Oc2cc(O)cc(O)c12)C')
    ]
    # Interflavan linkages patterns
    interflavan_b_type = Chem.MolFromSmarts('C1C(O)c2cc(O)ccc2OC1')
    interflavan_a_type = Chem.MolFromSmarts('OC1[C@H](OC2c3cc(O)ccc3C(OC2)C1)CO')

    # Check if at least one flavan-3-ol pattern matches
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavan3ol_patterns):
        return False, "No flavan-3-ol unit found"

    # Ensure either type of interflavan linkage is present
    if not (mol.HasSubstructMatch(interflavan_b_type) or mol.HasSubstructMatch(interflavan_a_type)):
        return False, "Insufficient evidence of distinctive proanthocyanidin linkages"

    # Check for presence of aromatic rings with at least two hydroxyls
    num_aromatic_hydroxyls = sum(mol.GetSubstructMatches(Chem.MolFromSmarts('c1c(O)cc(O)c(O)c1')))
    if num_aromatic_hydroxyls < 2:
        return False, "Aromatic hydroxyl content too low for proanthocyanidin"

    # Ensure chiral centers align with expected chemistry (at least 2)
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if num_chiral_centers < 2:
        return False, "Not enough chiral centers; proanthocyanidins typically have multiple"

    return True, "Contains features consistent with a proanthocyanidin: flavan-3-ol units and specific linkages"