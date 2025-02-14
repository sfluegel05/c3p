"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem
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

    # 1. Check for hydroperoxy group (-OOH or -OO-)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX1]")  # Matches -OOH and -OO-
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy group found"

    # 2. Check for a fatty acid chain (at least 12 carbons with a carboxyl group or OH)
    
    # Check for a carbonyl group (C=O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found, not a lipid."

    # Check for a carboxyl group or another OH
    carboxyl_or_oh_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    oh_pattern = Chem.MolFromSmarts("[OH1]")
    if not (mol.HasSubstructMatch(carboxyl_or_oh_pattern) or mol.HasSubstructMatch(oh_pattern)):
            return False, "No carboxyl group or OH found, not a lipid."


    # Check for a long carbon chain
    fatty_acid_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") # 10 C
    
    if not mol.HasSubstructMatch(fatty_acid_chain_pattern):
         return False, "Too short chain, not a fatty acid"

    # Check for a sufficient number of rotatable bonds to ensure long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too few rotatable bonds, not a long enough chain"

    # 3. Check for rings. Fatty acid derived hydroperoxides are typically acyclic
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings > 0:
        return False, "Contains rings, not a typical fatty acid hydroperoxide"

    return True, "Contains at least one hydroperoxy group and a fatty acid chain"