"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones typically have a benzoquinone core with methoxy groups and an isoprenyl side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzoquinone core with 2,3-dimethoxy groups
    benzoquinone_core_pattern = Chem.MolFromSmarts("COc1cc(=O)c(CC=C(C)C)cc1=O")
    if not mol.HasSubstructMatch(benzoquinone_core_pattern):
        return False, "No benzoquinone core with 2,3-dimethoxy groups found"

    # Assume long isoprenoid tail by checking for multiple repeating units R/C=C/R'
    isoprenoid_pattern = Chem.MolFromSmarts("C=C(C)C")
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_pattern)
    if len(isoprenoid_matches) < 2:
        return False, "Insufficient isoprenoid units"

    return True, "Contains benzoquinone core with methoxy groups and sufficient isoprenoid units"

# Test cases
test_smiles = [
    "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-9
    "O=C1C(OC)=C(OC)C(=O)C=C1C(C(O)C)C",  # 5-(3-hydroxybutan-2-yl)-2,3-dimethoxycyclohexa-2,5-diene-1,4-dione
]

for smiles in test_smiles:
    result, reason = is_ubiquinones(smiles)
    print(f"SMILES: {smiles}\nIs Ubiquinone: {result}\nReason: {reason}\n")