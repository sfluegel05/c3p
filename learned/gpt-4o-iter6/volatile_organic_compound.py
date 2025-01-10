"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is assumed to have a boiling point less than or equal to 250 degrees C, suggesting lower
    molecular weight, certain functional groups, and structural characteristics favoring volatility.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered a volatile organic compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count number of rotatable bonds
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Identify significant functional groups and simple structures
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    is_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    ester_pattern = Chem.MolFromSmarts("[$(OC(=O))$(C(=O)O)]")
    is_ester = mol.HasSubstructMatch(ester_pattern)
    simple_aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    is_simple_aromatic = mol.HasSubstructMatch(simple_aromatic_pattern)
    small_alkene_pattern = Chem.MolFromSmarts("C=C")
    is_small_alkene = mol.HasSubstructMatch(small_alkene_pattern)
    
    # Determine non-volatile functional groups and structural complexity
    high_boiling_pattern = Chem.MolFromSmarts("[$(NC=O),$(C(=O)N),$(S=O)]")
    has_high_boiling_features = mol.HasSubstructMatch(high_boiling_pattern)
    
    # Combine the rule-based approach
    if mol_wt <= 350 and c_count <= 20 and (is_alcohol or is_ester or is_small_alkene or is_simple_aromatic or n_rotatable_bonds <= 4) and not has_high_boiling_features:
        return True, f"Molecular weight {mol_wt}, carbon count {c_count}, functional groups, and structural properties suggest volatility."
    else:
        return False, f"Molecular weight {mol_wt}, carbon count {c_count}, and lack of indicative structural features contribute to non-volatility."