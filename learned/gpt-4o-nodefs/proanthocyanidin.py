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

    # Look for flavan-3-ol core pattern
    # Catechin and epicatechin are common units, requiring phenolic rings and multiple chiral centers
    flavan3ol_pattern = Chem.MolFromSmarts('Oc1cc(O)c2C[C@@H]([C@H](O)C2Oc2cc(O)cc(O)c12)C')
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "No flavan-3-ol unit found"

    # Check for interflavan linkages patterns, including both link types A and B
    # B-type C4-C8 linkage, A-type involves C4-C8 and additional ether linkages
    interflavan_b_type = Chem.MolFromSmarts('C1C(O)c2cc(O)ccc2OC1')
    interflavan_a_type = Chem.MolFromSmarts('OC1[C@H](OC2c3cc(O)ccc3C(OC2)C1)CO')

    if not (mol.HasSubstructMatch(interflavan_b_type) or mol.HasSubstructMatch(interflavan_a_type)):
        return False, "Insufficient evidence of distinctive proanthocyanidin linkages"

    # Ensure aromatic rings with hydroxyl groups exist
    num_aromatic_oh = len([atom for atom in mol.GetAtoms() 
                           if atom.GetIsAromatic() and atom.GetSymbol() == 'O'])
    if num_aromatic_oh < 2:
        return False, "Aromatic hydroxyl content too low for proanthocyanidin"

    # Verify the presence of chiral centers
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if num_chiral_centers < 2:
        return False, "Not enough chiral centers; proanthocyanidins typically have multiple"

    return True, "Contains features consistent with a proanthocyanidin: flavan-3-ol units and specific linkages"