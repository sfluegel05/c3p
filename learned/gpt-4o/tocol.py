"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton that is substituted at position 2 
    by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern for chromanol core: aromatic benzopyran with hydroxyl group
    # This pattern captures chroman-6-ol skeleton
    chromanol_pattern = Chem.MolFromSmarts("Oc1ccccc1C2CCCCO2")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chromanol core found"
    
    # Look for substituent at position 2 with isoprenoid characteristics.
    # Assumes the side chain is a long alkyl/isoprenoid without relying on specific stereochemistry
    hydrocarbon_substituent_pattern = Chem.MolFromSmarts("C(C)(C)C")
    isoprenoid_matches = mol.GetSubstructMatches(hydrocarbon_substituent_pattern)
    if len(isoprenoid_matches) < 3:
        return False, "Insufficient isoprenoid features in substituent"

    return True, "Contains chromanol core and appropriate hydrocarbon substitution with isoprenoid units"