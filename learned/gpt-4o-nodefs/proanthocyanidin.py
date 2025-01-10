"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    Proanthocyanidins are polyphenolic structures characterized by flavan-3-ol units
    linked together, featuring aromatic hydroxyl groups and structures like catechin/epicatechin.

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

    # Look for flavan-3-ol core pattern
    flavan3ol_pattern = Chem.MolFromSmarts('OC1C(COC2=CC(O)=C(O)C=C2)C3=CC(O)=C(O)C=C13')
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "No flavan-3-ol unit found"

    # Check for linkage between flavan-3-ol units (common proanthocyanidin linkages)
    # This pattern checks for interflavan linkages like B-type linkages
    interflavan_pattern = Chem.MolFromSmarts('OC1C(O)C2=CC(O)=CC(O)=C2C1C')
    if not mol.HasSubstructMatch(interflavan_pattern):
        return False, "Insufficient evidence of proanthocyanidin linkages"

    # Ensure multiple chiral centers are present; proanthocyanidins contain them.
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if num_chiral_centers < 2:
        return False, "Not enough chiral centers; proanthocyanidins typically have multiple"

    return True, "Contains features consistent with a proanthocyanidin: flavan-3-ol units linked with specific proanthocyanidin linkages"