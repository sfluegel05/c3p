"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

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

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Pattern for a 12-carbon saturated chain attached to an ester
    # The pattern ensures we find exactly 12 carbons in a row
    laurate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]")
    
    if mol.HasSubstructMatch(laurate_pattern):
        # For each match, verify it's a proper laurate ester
        matches = mol.GetSubstructMatches(laurate_pattern)
        for match in matches:
            # Get the carbons in the chain
            chain_carbons = match[2:]  # Skip O and C(=O)
            
            # Verify these carbons form a linear chain
            is_linear = True
            for i, c_idx in enumerate(chain_carbons):
                atom = mol.GetAtomWithIdx(c_idx)
                # Check correct number of connections
                if i == 0:  # First carbon
                    if len([n for n in atom.GetNeighbors() if n.GetIdx() not in match]) > 0:
                        is_linear = False
                        break
                elif i == len(chain_carbons) - 1:  # Last carbon (methyl)
                    if atom.GetDegree() != 1:
                        is_linear = False
                        break
                else:  # Middle carbons
                    if len([n for n in atom.GetNeighbors() if n.GetIdx() not in chain_carbons]) > 0:
                        is_linear = False
                        break
                
                # Verify no double bonds in the chain
                for bond in atom.GetBonds():
                    if bond.GetBondType() != Chem.BondType.SINGLE:
                        is_linear = False
                        break
            
            if is_linear:
                return True, "Contains dodecanoate (laurate) ester group with linear saturated C12 chain"
    
    return False, "No dodecanoate ester group with appropriate chain length found"