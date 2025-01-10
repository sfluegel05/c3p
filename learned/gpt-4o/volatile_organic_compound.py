"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is likely to be a volatile organic compound (VOC)
    based on its SMILES string. This is a heuristic guess considering molecular
    weight, composition, and presence of volatile functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is possibly a VOC, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Determine molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt > 250:
        return False, f"Molecular weight too high ({mol_wt:.2f} Da) for typical VOC"

    # Check for non-metallic typical VOC composition
    non_metallic = all(atom.GetAtomicNum() in {5, 6, 7, 8, 9, 15, 16, 17, 35, 53} for atom in mol.GetAtoms())
    if not non_metallic:
        return False, "Contains elements not typical of VOCs"

    # Identify presence of volatile functional groups: alcohols, ethers, aldehydes, ketones
    volatile_functional_groups = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Aldehyde
        Chem.MolFromSmarts("[#6]-[#8]-[#6]"),    # Ether
        Chem.MolFromSmarts("[CX3](=O)[#6]"),     # Ketone
        Chem.MolFromSmarts("[#6][CX2]([#6])([#6])[#8X2H]"),  # Alcohol
    ]
    for pattern in volatile_functional_groups:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains typical volatile functional groups"

    # Calculate descriptors and other checks (additional heuristics)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable_bonds > 6:
        return False, "Too many rotatable bonds for a VOC"

    return True, "Likely a VOC based on molecular weight, typical elements, and structure"

# __metadata__ or additional configuration can be added here if needed.