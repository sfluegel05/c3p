"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: CHEBI:28425 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is a phosphatidylglycerol composed of two phosphatidic acid molecules linked to a central glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly two phosphorus atoms (from two phosphate groups)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 2:
        return False, f"Found {p_count} phosphorus atoms, need exactly 2"

    # Corrected ester group detection (R-O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)[#6]")  # Matches O-C(=O)-C
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 4"

    # Generalized glycerol backbone pattern (3 carbons with oxygen substituents)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][O]-[CH][O]-[CH2][O]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) < 3:
        return False, f"Found {len(glycerol_matches)} glycerol backbones, need at least 3"

    # Check phosphate group presence
    phosphate_pattern = Chem.MolFromSmarts("[O]P(=O)([O])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Insufficient phosphate groups"

    # Check molecular weight (cardiolipins are typically >1200 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 1200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    # Check for long carbon chains (at least 12 carbons in each fatty acid)
    chain_pattern = Chem.MolFromSmarts("[#6]C(=O)O[#6]")  # Carbon adjacent to ester carbonyl
    for match in mol.GetSubstructMatches(chain_pattern):
        chain_start = match[0]
        atom = mol.GetAtomWithIdx(chain_start)
        chain_length = 0
        
        # Traverse the chain backwards (away from carbonyl)
        while True:
            neighbors = [n for n in atom.GetNeighbors() 
                        if n.GetAtomicNum() == 6 
                        and n.GetIdx() != match[1]]  # Avoid going back to carbonyl
            if len(neighbors) != 1:
                break
            chain_length += 1
            atom = neighbors[0]
            
        if chain_length < 12:
            return False, f"Fatty acid chain too short ({chain_length} carbons)"

    return True, "Contains two phosphatidic acids linked via central glycerol with two phosphate groups"