"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are derived from prostanoic acid, a C20 fatty acid with a cyclopentane ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a substituted cyclopentane or related ring
    # Allow for 5-membered rings with C and/or O, double bonds, and substitutions
    cyclopentane_pattern = Chem.MolFromSmarts("[C,O]1([C,O])([C,O])([C,O])([C,O])1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
       return False, "No cyclopentane or related ring found"

    # 2. Check for a carboxyl acid, ester, or amide group
    acid_ester_amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])O") # -C(=O)O, or -C(=O)OR
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    if not (mol.HasSubstructMatch(acid_ester_amide_pattern) or mol.HasSubstructMatch(amide_pattern)):
        return False, "No carboxyl, ester, or amide group found"
        
    #3. Check for presence of carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No carbon chain found"

    #4. Check for at least 2 alcohol groups
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]") # -OH
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) < 2:
        return False, f"Must have at least 2 alcohol groups"

    #5. Check for at least 3 oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Must have at least 3 oxygens, found {o_count}"

    # 6. Check for a carbon chain connected to the cyclopentane ring and a carbonyl group
    # This should help prevent false positives. This check is similar to the core pattern check from before
    core_pattern = Chem.MolFromSmarts("[C,O]1([C,O])([C,O])([C,O])([C,O])1~[C]~[CX3]=[OX1]")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core prostaglandin pattern not found"

    return True, "Meets basic prostaglandin criteria"