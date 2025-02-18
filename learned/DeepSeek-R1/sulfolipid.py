"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:73404 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid group (C-SO3) attached via a carbon-sulfur bond to a lipid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfonic acid group (C-SO3H)
    sulfonic_acid_pattern = Chem.MolFromSmarts("[C][S](=O)(=O)[O]")
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group (C-SO3H) found"

    # Check lipid characteristics through carbon count and molecular weight
    # Count total carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Insufficient carbons ({c_count}) for lipid structure"

    # Verify significant molecular weight (>500 Da typical)
    mol_wt = AllChem.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for lipid"

    return True, "Contains sulfonic acid group attached to lipid structure"