"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Chromanol core pattern SMARTS: includes benzopyran-6-ol structure   
    chromanol_pattern = Chem.MolFromSmarts("Oc1ccc(CC2CCCCO2)c(c)c1")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chromanol core found"

    # Isoprenoid chain characteristics using common isoprene unit structure (C5H8)
    # Yet flexible enough to allow saturation and unsaturation as permitted
    isoprenoid_unit = Chem.MolFromSmarts("C(C)(C)C")
    num_isoprenoid_units = len(mol.GetSubstructMatches(isoprenoid_unit))
    
    if num_isoprenoid_units < 3:
        return False, f"Insufficient isoprenoid features: {num_isoprenoid_units} found, need at least 3"

    return True, "Contains chromanol core and appropriate hydrocarbon substitution with isoprenoid units"