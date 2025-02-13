"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Classifies a SMILES string as a triglyceride if it fulfills the structural requirements:
    1. A glycerol backbone.
    2. Three ester groups resulting from the condensation with fatty acids.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if it is a triglyceride, False otherwise
        str: Explanation of the result
    """
    
    # Parse SMILES to obtain molecular representation
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect glycerol-like backbone: C-C-C 3 oxygens connected (indicative of glycerol)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Detect 3 ester linkages: -O-C(=O)-
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_count = len(mol.GetSubstructMatches(ester_pattern))
    if ester_count != 3:
        return False, f"Found {ester_count} ester groups, expected 3"

    # Calculate number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, "Insufficient carbon count for typical triglycerides"

    # Calculate molecular weight, typical triglyceride is >500 g/mol
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 500:
        return False, "Molecular weight too low for triglyceride"

    return True, "SMILES represents a triglyceride with a glycerol backbone and three ester groups"