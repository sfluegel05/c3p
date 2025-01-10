"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC typically has an initial boiling point less than or equal to 250 degrees C at 101.3 kPa,
    characterized by specific structural characteristics favoring volatility.

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
    
    # Consider known volatile structures like cyclic ethers or small lactones
    small_lactone_pattern = Chem.MolFromSmarts("C1OC(=O)CC1")
    is_small_lactone = mol.HasSubstructMatch(small_lactone_pattern)
    
    small_cyclic_alcohol_pattern = Chem.MolFromSmarts("C1(CCCCC1)O")
    is_small_cyclic_alcohol = mol.HasSubstructMatch(small_cyclic_alcohol_pattern)

    # Updated rules for classification based on refined criteria
    if (((mol_wt <= 300 and c_count <= 20) or (is_alcohol and c_count <= 15) or 
         (is_alkene and c_count <= 15) or (is_halogenated and c_count <= 12) or 
         is_small_lactone or is_small_cyclic_alcohol) and 
        n_rotatable_bonds <= 6):
        
        return True, f"Molecular weight {mol_wt}, carbon count {c_count}, functional groups, and structural properties suggest volatility."
    
    return False, f"Molecular weight {mol_wt}, carbon count {c_count}, and structural characteristics suggest non-volatility."