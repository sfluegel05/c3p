"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:16676 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.
    The term was originally restricted to natural products, but is now also used to describe
    semi-synthetic and fully synthetic compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 1-benzopyran core or related fused ring systems
    benzopyran_pattern = Chem.MolFromSmarts("O=C1C=2C=CC=CC2=CC=C1")
    benzopyran_matches = mol.GetSubstructMatches(benzopyran_pattern)
    if not benzopyran_matches:
        return False, "No 1-benzopyran core or related fused ring system found"

    # Look for aryl substituent at position 3 (allowing for substituents)
    aryl_pattern = Chem.MolFromSmarts("[a;!r3]")
    aryl_matches = mol.GetSubstructMatches(aryl_pattern, maxMatches=1, useChirality=True)
    if not aryl_matches:
        return False, "No aryl substituent found at position 3"
    
    aryl_idx = aryl_matches[0][0]
    atom = mol.GetAtomWithIdx(aryl_idx)
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    if 2 not in neighbors:
        return False, "Aryl substituent not attached at position 3"

    # Check for additional substituents on the benzopyran core or aryl ring
    allowed_substituents = ["C", "O", "N", "S", "Cl", "Br", "F", "I"]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [6, 7, 8, 9, 16, 17, 35, 53] and atom.GetTotalNumHs() > 0:
            if atom.GetSymbol() not in allowed_substituents:
                return False, f"Unexpected substituent: {atom.GetSymbol()}"

    # Check for typical molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical range for isoflavonoids"

    return True, "Meets structural criteria for isoflavonoids"