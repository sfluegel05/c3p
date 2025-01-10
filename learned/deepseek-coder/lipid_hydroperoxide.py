"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: CHEBI:36044 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved hydroperoxy group pattern (-OOH)
    # This pattern matches hydroperoxy groups attached to any carbon
    hydroperoxy_pattern = Chem.MolFromSmarts("[CX4][OX2][OX1]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    
    if not hydroperoxy_matches:
        # Try alternative pattern for hydroperoxy groups attached to double-bonded carbons
        hydroperoxy_pattern = Chem.MolFromSmarts("[CX3]=[CX3][OX2][OX1]")
        hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
        if not hydroperoxy_matches:
            return False, "No hydroperoxy group found"

    # Basic lipid check: at least 10 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons for lipid"

    # Check molecular weight - lipids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for lipid"

    return True, "Contains hydroperoxy group and meets lipid criteria"