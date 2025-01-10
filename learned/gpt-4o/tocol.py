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

    # Chromanol 6-ol structure is characterized typically by a 2H-1-benzopyran-6-ol structure
    chromanol_core = Chem.MolFromSmarts("Oc1ccc2OCC(Cc2c1)C")
    if not mol.HasSubstructMatch(chromanol_core):
        return False, "No chromanol core (2H-benzopyran-6-ol) found"

    # Isoprenoid chains are built from repeated C5 units.
    # Here, we're simplifying by checking for broad isoprenoid patterns.
    isoprenoid_unit = Chem.MolFromSmarts("C(C)(C)CC")
    num_isoprenoid_units = len(mol.GetSubstructMatches(isoprenoid_unit))
    
    # Assume that variations in the chain (i.e., degree of unsaturation) mean we should at least look for three chains, base
    # This logic is flexible and some manual tuning might be needed for edge cases not covered.
    if num_isoprenoid_units < 3 and not (mol.HasSubstructMatch(Chem.MolFromSmarts("C=C"))):
        return False, f"Insufficient isoprenoid units: {num_isoprenoid_units} found, with minimal unsaturation"

    return True, "Contains chromanol core with appropriate hydrocarbon chain substitution"