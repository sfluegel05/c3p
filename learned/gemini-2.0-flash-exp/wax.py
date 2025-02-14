"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are generally long-chain esters or alcohols with flexibility
    and moderate to high molecular weight.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for long carbon chains (>= 20 C)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Too few carbons ({c_count}), typical wax should have at least 20."

    # 2. Check for ester linkages (or long-chain alcohols).
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H1]") #check for long chain alcohols
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    if len(ester_matches) == 0 and len(alcohol_matches) == 0:
        return False, f"No ester groups or long-chain alcohol detected. Wax should have one of these."

    #3 Check number of oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if len(ester_matches) > 0 and o_count < len(ester_matches)*2 : #if there is an ester, must have > 2 oxygen
        return False, "Inconsistent number of oxygens for number of esters"
    if len(alcohol_matches) > 0 and o_count < len(alcohol_matches) : #if there is an alcohol, must have > 1 oxygen
        return False, "Inconsistent number of oxygens for number of alcohols"

    # 4. Check for rotatable bonds (at least 8-10)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, f"Too few rotatable bonds ({n_rotatable}). Should be at least 8"

    # 5. Check molecular weight (approx >= 300)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt}). Should be at least 300"

    # 6. Check number of double bonds
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) > 0:
        pass #double bonds are good
    
    return True, "Meets the criteria for a wax (long chains, ester or alcohol, rotatable bonds, good MW)."