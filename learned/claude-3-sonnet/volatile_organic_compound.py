"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:51564 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from typing import Tuple

def is_volatile_organic_compound(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is defined as an organic compound with an initial boiling point <= 250°C at 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule is organic (contains carbon)
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
        return False, "Molecule does not contain carbon, not an organic compound"
    
    # Look for functional groups and substructures commonly found in VOCs
    voc_patterns = (
        "[OX1]",                   # Aldehydes and ketones
        "[NX1]#[CX2]",             # Nitriles
        "[CX3]=[CX3]",             # Alkenes
        "[OX2][CX3]=[OX1]",        # Esters and carboxylic acids
        "[CX4]=[OX1]",             # Aldehydes
        "[NX3][CX3]=[OX1]",        # Amides
        "[SX2]=[OX1]",             # Sulfoxides and sulfones
        "[CX4][OX2][CX3]"          # Alcohols
    )
    
    for pattern in voc_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            # Estimate boiling point using a combination of descriptors
            mol_wt = Descriptors.ExactMolWt(mol)
            n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            est_boiling_point = 198.18 + 20.07 * mol_wt - 7.64 * n_rotatable - 6.8 * tpsa
            
            if est_boiling_point <= 250:
                return True, f"Estimated boiling point of {est_boiling_point:.2f}°C is <= 250°C"
            else:
                return False, f"Estimated boiling point of {est_boiling_point:.2f}°C is > 250°C"
    
    return False, "No functional groups or substructures characteristic of VOCs found"