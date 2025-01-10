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
        bool: True if the molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS pattern for flavan-3-ol backbone (C6-C3-C6 with a hydroxyl on C3)
    flavan_3_ol_pattern = Chem.MolFromSmarts("C1C(C2=C(O)C=CC(O)=C2)O[C@H]([C@@H]1)c3ccc(O)cc3")
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol backbone found"
    
    # Counting hydroxyl and other typical substituents (flexible substitution pattern)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    methoxy_pattern = Chem.MolFromSmarts("CO")
    ester_pattern = Chem.MolFromSmarts("O=C(O)")
    
    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    num_oxygen_containing_subs = len(hydroxyl_matches) + len(methoxy_matches) + len(ester_matches)
    
    if num_oxygen_containing_subs < 2:
        return False, f"Insufficient oxygen-containing groups, found {num_oxygen_containing_subs}"

    # Verify stereochemistry if there are chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) > 0 and len(chiral_centers) < 2:
        return False, f"Insufficient chiral centers, found {len(chiral_centers)}"

    # Assume various substitutions while maintaining flavan-3-ol identity
    return True, "Contains flavan-3-ol backbone with sufficient oxygen-containing groups and stereochemistry"