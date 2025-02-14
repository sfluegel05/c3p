"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is a fatty acid that contains one or more rings as part of the fatty acid chain.

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

    # Check for carboxylic acid or ester group directly connected to a carbon chain
    acid_or_ester_pattern = Chem.MolFromSmarts('[CX3,CX4](=[OX1])[OX2]')
    acid_or_ester_matches = mol.GetSubstructMatches(acid_or_ester_pattern)

    if not acid_or_ester_matches:
        return False, "No carboxylic acid or ester group found"


    # Check for a carbon chain, containing a ring, connected to acid group
    # Modified to look for chain with at least 3 carbons including one carbon that is part of a ring
    # with the carboxyl group attached to a carbon of the chain.
    # This pattern looks for the carboxyl carbon atom, then a carbon connected to it, then a carbon inside a ring and then a carbon atom connected to it,
    # to ensure that the ring is part of the chain.
    
    cyclic_chain_pattern = Chem.MolFromSmarts('[CX3,CX4](=[OX1])[OX2][#6]-[#6;R]-[#6]')

    
    found_chain_match = False
    for match in acid_or_ester_matches:
        for atom_index in match: #check if the carbon of the acid/ester is part of a chain
            if mol.GetAtomWithIdx(atom_index).GetAtomicNum() == 6:
                for neighbor in mol.GetAtomWithIdx(atom_index).GetNeighbors(): #find carbons attached to carboxyl
                    if neighbor.GetAtomicNum() == 6:
                        #check for full chain pattern with carboxyl group:
                       matches = mol.GetSubstructMatches(cyclic_chain_pattern)
                       if matches:
                        for m in matches:
                            if atom_index in m:
                             found_chain_match = True
                               

    if not found_chain_match:
      return False, "No ring found as part of a fatty acid carbon chain"

    #count carbons in whole molecule to ensure it's at least a short fatty acid (5 or more)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbon atoms for a fatty acid"

    return True, "Contains ring and fatty acid substructures"