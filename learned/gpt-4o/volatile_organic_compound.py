"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is likely to be a volatile organic compound (VOC)
    based on its SMILES string. Since boiling point cannot be directly determined,
    this function is a heuristic guess based on molecular weight and composition.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is possibly a VOC, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use molecular weight as a proxy: most VOCs are lighter
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt > 300:
        return False, f"Molecular weight too high ({mol_wt:.2f} Da) for typical VOC"

    # Check molecular composition: presence of non-metal elements typical for organics
    non_metallic = all(atom.GetAtomicNum() in {5, 6, 7, 8, 9, 15, 16, 17, 35, 53} for atom in mol.GetAtoms())
    if not non_metallic:
        return False, "Contains elements not typical of VOCs"
    
    return True, "Likely a VOC based on molecular weight and typical elements"

# __metadata__ or additional configuration can be added here if needed.