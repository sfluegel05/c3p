"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: CHEBI:17411 lysophosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are monoacylglycerol phosphates with one acyl chain and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one ester group (acyl chain)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3]=[OX1]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for phosphate group connected to carbon backbone
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)(O)(O)O[C]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group connected to carbon backbone"

    # Get ester and phosphate attachment points
    try:
        # Ester oxygen and carbonyl carbon
        ester_o_idx = ester_matches[0][0]
        ester_c_idx = ester_matches[0][1]
        
        # Phosphate oxygen and attached carbon (glycerol backbone)
        phosphate_p_idx = phosphate_matches[0][0]
        phosphate_o_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(phosphate_p_idx).GetNeighbors() 
                          if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(phosphate_p_idx, n.GetIdx()).GetBondType() == Chem.BondType.SINGLE][0]
        phosphate_c_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(phosphate_o_idx).GetNeighbors() if n.GetAtomicNum() == 6][0]

        # Check if ester and phosphate are on a three-carbon backbone
        path = Chem.GetShortestPath(mol, ester_c_idx, phosphate_c_idx)
        if len(path) not in [2, 3]:  # Allow adjacent or separated by one carbon
            return False, "Ester and phosphate not on three-carbon backbone"
            
        # Verify all atoms in path are carbons
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path):
            return False, "Backbone contains non-carbon atoms"

    except:
        return False, "Error analyzing molecular structure"

    # Optional: Check acyl chain length (at least 8 carbons)
    try:
        acyl_chain = Chem.MolFromSmarts(f"[#{ester_c_idx}]-[OX2]-[CX3]=O")
        chain_length = len(Chem.MolFromSmarts("[C]").GetSubstructMatches(mol.GetAtomWithIdx(ester_c_idx).GetNeighbors()[0].GetOwningMol()))
        if chain_length < 8:
            return False, f"Acyl chain too short ({chain_length} carbons)"
    except:
        pass  # Skip chain length check if analysis fails

    return True, "Monoacylglycerol phosphate with one acyl group and one phosphate group on three-carbon backbone"