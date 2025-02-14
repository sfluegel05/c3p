"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is a fatty acid that contains one or more rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ring structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No ring structure found"

    # Check for carboxylic acid or ester group directly connected to a carbon chain, and a long carbon chain using SMARTS
    # Carboxylic acid pattern and ester group are allowed
    acid_or_ester_pattern = Chem.MolFromSmarts('[CX4,CX3](=[OX1])([OX2])[#6]') # carbonyl attached to O, attached to carbon
    acid_or_ester_matches = mol.GetSubstructMatches(acid_or_ester_pattern)

    if not acid_or_ester_matches:
         return False, "No carboxylic acid or ester group found attached to a carbon chain"

    # Check for a hydrocarbon chain attached to the carboxylic acid group with at least 4 carbons
    hydrocarbon_chain_pattern = Chem.MolFromSmarts('[#6][#6]([#6])[#6]')
    
    found_long_chain = False
    for match in acid_or_ester_matches:
        for atom_index in match:
             atom = mol.GetAtomWithIdx(atom_index)
             if atom.GetAtomicNum() == 6:
                #check if this carbon is part of the hydrocarbon chain by looking at its neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        for neighbor2 in neighbor.GetNeighbors():
                            if neighbor2.GetAtomicNum() == 6 and neighbor2.GetIdx() != atom_index :
                                for neighbor3 in neighbor2.GetNeighbors():
                                    if neighbor3.GetAtomicNum() == 6 and neighbor3.GetIdx() != neighbor.GetIdx():
                                        found_long_chain=True
    if not found_long_chain:
      return False, "No long hydrocarbon chain attached to the acid/ester"

    # Check if the ring is part of the fatty acid chain (not a side group)
    # Here we're checking if a carbon from the ring is directly attached to a carbon chain and the acid group
    found_ring_chain = False
    for ring_atom_idx in ring_info.AtomRings():
        for atom_idx in ring_atom_idx:
             atom = mol.GetAtomWithIdx(atom_idx)
             if atom.GetAtomicNum() == 6:
                # Check if this carbon is attached to a carbon of a long chain
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        for neighbor2 in neighbor.GetNeighbors():
                            if neighbor2.GetAtomicNum() == 6 and neighbor2.GetIdx() != atom_idx :
                                found_ring_chain = True
    if not found_ring_chain:
         return False, "Ring is not part of the fatty acid chain."

    # count carbons: fatty acids typically have >4 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbon atoms for a fatty acid"
   
    return True, "Contains ring and fatty acid substructures"