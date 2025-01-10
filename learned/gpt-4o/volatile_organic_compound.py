"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is likely to be a volatile organic compound (VOC)
    based on its SMILES string, considering molecular weight, presence of volatile
    functional groups, and other chemical features.

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
    if mol_wt > 300:
        return False, f"Molecular weight too high ({mol_wt:.2f} Da) for typical VOC"

    # Check for elements common in VOCs including halogenated ones
    non_metallic = all(atom.GetAtomicNum() in {1, 5, 6, 7, 8, 9, 17, 35, 53} for atom in mol.GetAtoms())
    if not non_metallic:
        return False, "Contains elements not typical of VOCs"
    
    # Identify presence of volatile functional groups: alcohols, ethers, aldehydes, ketones, chloroalkanes
    volatile_functional_groups = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Aldehydes
        Chem.MolFromSmarts("[#6]-[#8]-[#6]"),    # Ethers
        Chem.MolFromSmarts("[CX3](=O)[#6]"),     # Ketones
        Chem.MolFromSmarts("[#6][CX2]([#6])([#6])[#8X2H]"),  # Alcohols
        Chem.MolFromSmarts("[Cl,Br,F,I][#6]"),   # Halogenated hydrocarbons
    ]
    for pattern in volatile_functional_groups:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains typical volatile functional groups"
    
    # Calculate number of rotatable bonds
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)

    # Apply a rotable bond and ring heuristic (favor cyclic structures)
    if num_rotatable_bonds + num_rings > 8:
        return False, "Too many unbounded or undefined bonds for a VOC"

    return True, "Likely a VOC based on molecular weight, typical elements, and structure"

# __metadata__ or additional configuration can be added here if needed.