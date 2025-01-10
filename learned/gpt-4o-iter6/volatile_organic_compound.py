"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is assumed to commonly have a boiling point less than or equal to 250 degrees C,
    suggesting they often have specific structural characteristics favoring volatility.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is considered a VOC, False otherwise.
        str: Reason for classification.
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
    
    # Identify functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    is_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    alkene_pattern = Chem.MolFromSmarts("C=C")
    is_alkene = mol.HasSubstructMatch(alkene_pattern)
    halogenated_pattern = Chem.MolFromSmarts("[CX4][F,Cl,Br,I]")
    is_halogenated = mol.HasSubstructMatch(halogenated_pattern)
    simple_alkane_pattern = Chem.MolFromSmarts("[CX4;!R][CX4;!R]")
    is_simple_alkane = mol.HasSubstructMatch(simple_alkane_pattern)
    
    # Rules for classification based on refined criteria
    if (mol_wt <= 350 and c_count <= 20 and
        (is_alcohol or is_alkene or is_simple_alkane or is_halogenated or n_rotatable_bonds <= 7)):
        return True, f"Molecular weight {mol_wt}, carbon count {c_count}, functional groups, and structural properties suggest volatility."
    
    return False, f"Molecular weight {mol_wt}, carbon count {c_count}, and structural characteristics suggest non-volatility."