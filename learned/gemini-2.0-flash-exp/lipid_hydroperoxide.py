"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is a lipid containing one or more hydroperoxy (-OOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for hydroperoxy group (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX1]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"

    # 2. Check for a long carbon chain (at least 4 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
         return False, "No long carbon chain"
    
    # Check for rotatable bonds (at least 3 to imply the chain)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Chain too short (few rotatable bonds)."

    # 3. (Optional) Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200: # common lipids typically have MW higher than 200. Can be lowered if needed
        return False, "Molecular weight too low for a lipid"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons for lipid"

    return True, "Contains a hydroperoxy group and long carbon chain"