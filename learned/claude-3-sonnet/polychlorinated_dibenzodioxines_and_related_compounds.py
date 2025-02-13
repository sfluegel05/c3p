"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: CHEBI:51881 polychlorinated dibenzodioxins and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin or a related compound
    based on its SMILES string. Related compounds include polychlorinated dibenzofurans,
    polychlorinated biphenyls, and polybrominated biphenyls.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorinated dibenzodioxin or related compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for polychlorinated dibenzodioxins and related compounds
    pcdibenzodioxins_pattern = Chem.MolFromSmarts("*c1c(Cl)c(Cl)c2Oc3c(Cl)c(Cl)c(Cl)c(Cl)c3Oc2c1*")
    pcdibenzofurans_pattern = Chem.MolFromSmarts("*c1c(Cl)c(Cl)c2oc3c(Cl)c(Cl)c(Cl)c(Cl)c3c2c1*")
    pcbiphenyls_pattern = Chem.MolFromSmarts("Clc1ccc(c*c1Cl)-c1c(Cl)ccc(Cl)c1Cl")
    pbbiphenyls_pattern = Chem.MolFromSmarts("Brc1ccc(c*c1Br)-c1c(Br)ccc(Br)c1Br")

    # Check for polychlorinated dibenzodioxins
    if mol.HasSubstructMatch(pcdibenzodioxins_pattern):
        return True, "Polychlorinated dibenzodioxin"

    # Check for polychlorinated dibenzofurans
    if mol.HasSubstructMatch(pcdibenzofurans_pattern):
        return True, "Polychlorinated dibenzofuran"

    # Check for polychlorinated biphenyls
    if mol.HasSubstructMatch(pcbiphenyls_pattern):
        return True, "Polychlorinated biphenyl"

    # Check for polybrominated biphenyls
    if mol.HasSubstructMatch(pbbiphenyls_pattern):
        return True, "Polybrominated biphenyl"

    # If none of the patterns match, it's not a polychlorinated dibenzodioxin or related compound
    return False, "Not a polychlorinated dibenzodioxin or related compound"