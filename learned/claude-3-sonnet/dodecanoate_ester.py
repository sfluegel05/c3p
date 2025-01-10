"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester contains a lauric acid (C12:0) component attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Look for C12 saturated chain attached to ester
    # Pattern: ester carbon followed by 11 more carbons in a chain
    dodecanoyl_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]")
    
    if not mol.HasSubstructMatch(dodecanoyl_pattern):
        return False, "No dodecanoyl (C12:0) chain found attached to ester"
    
    # Get all matches of the dodecanoyl pattern
    matches = mol.GetSubstructMatches(dodecanoyl_pattern)
    
    # For each match, verify it's a saturated chain by checking no double bonds
    for match in matches:
        # Get the carbons in the chain (excluding the ester carbonyl)
        chain_atoms = match[2:]  # Skip O and C(=O)
        
        # Check for any double bonds in the chain
        has_double_bond = False
        for atom_idx in chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    has_double_bond = True
                    break
        
        # If we found a saturated chain, it's a dodecanoate ester
        if not has_double_bond:
            return True, "Contains dodecanoate (laurate) ester group"
    
    return False, "No saturated dodecanoyl chain found"