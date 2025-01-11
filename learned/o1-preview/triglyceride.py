"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is a glycerol backbone with three fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with three oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO*)(CO*)CO*")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for three ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 3"

    # Check for long carbon chains (fatty acid chains)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]COC(=O)[C@H]CCCC")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, "Less than 3 fatty acid chains attached via ester bonds"

    # Check molecular weight - triglycerides typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a triglyceride"

    return True, "Contains glycerol backbone with 3 fatty acid chains attached via ester bonds"