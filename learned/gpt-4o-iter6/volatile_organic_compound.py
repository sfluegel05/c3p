"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is assumed to have a boiling point less than or equal to 250 degree C, suggesting lower molecular weight
    and structural characteristics that favor volatility.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered a volatile organic compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Count number of rotatable bonds
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    # Determine presence of functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    is_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    ester_pattern = Chem.MolFromSmarts("[$(OC(=O))$(C(=O)O)]")
    is_ester = mol.HasSubstructMatch(ester_pattern)
    simple_aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    is_simple_aromatic = mol.HasSubstructMatch(simple_aromatic_pattern)
    
    # Rule-based combination
    if mol_wt <= 300 and (
        c_count <= 20 or (
        is_alcohol or is_ester or is_simple_aromatic or
        n_rotatable_bonds > 1
    )):
        return True, f"Molecular weight {mol_wt}, carbon count {c_count}, and structural features suggest potential volatility."
    else:
        return False, f"Molecular weight {mol_wt}, carbon count {c_count}, and lack of indicative structural features suggest non-volatility."

# Test the function with a SMILES string
# Example: 1-dodecene (SMILES: CCCCCCCCCCC=C), expected to be a VOC
example_smiles = "CCCCCCCCCCC=C"
is_volatile, reason = is_volatile_organic_compound(example_smiles)
print(is_volatile, reason)