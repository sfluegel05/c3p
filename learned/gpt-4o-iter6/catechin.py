"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are characterized by a flavan-3-ol backbone and substituted hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Catechin is a hydroxyflavan: Check for flavan-3-ol structure
    flavan_3_ol_pattern = Chem.MolFromSmarts("C1[C@@H](OC2=CC=C(O)C=C2)C(=O)CC3=CC(O)=CC(O)=C13")
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol backbone found"
    
    # Check for hydroxyl substitutions
    hydroxy_pattern = Chem.MolFromSmarts("c(O)c")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxyl_matches)}"

    # Verify stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 2:
        return False, f"Insufficient chiral centers, found {len(chiral_centers)}"

    return True, "Contains flavan-3-ol backbone with appropriate hydroxyl groups and stereochemistry"