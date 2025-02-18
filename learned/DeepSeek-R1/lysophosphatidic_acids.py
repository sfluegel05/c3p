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
    Determines if a molecule is a lysophosphatidic acid (LPA) based on its SMILES string.
    LPA must have:
    - Glycerol backbone (3 carbons)
    - Exactly one ester-linked acyl chain
    - Phosphate group on the backbone
    - Both groups attached to distinct ends of the glycerol chain
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Find ester groups (O-C=O)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]C(=O)"))
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups (need 1)"

    # Find phosphate groups (must be O-P=O connected to carbon chain)
    phosphate_matches = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[PX4](=O)(O)(O)OC")
    )
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Get attachment points
    try:
        # Ester: oxygen connected to carbonyl carbon
        ester_o = ester_matches[0][0]
        ester_c = ester_matches[0][1]

        # Phosphate: oxygen connected to carbon backbone
        phosphate_p = phosphate_matches[0][0]
        phosphate_o = [n.GetIdx() for n in mol.GetAtomWithIdx(phosphate_p).GetNeighbors()
                      if n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 0][0]
        phosphate_c = [n.GetIdx() for n in mol.GetAtomWithIdx(phosphate_o).GetNeighbors()
                      if n.GetAtomicNum() == 6][0]

        # Check backbone connection
        path = Chem.GetShortestPath(mol, ester_c, phosphate_c)
        if len(path) != 3:  # Should form C1-C2-C3 pattern
            return False, f"Ester/phosphate not on 3-carbon backbone (path length {len(path)})"

        # Verify all backbone atoms are carbons
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in path):
            return False, "Backbone contains non-carbon atoms"

        # Check for remaining hydroxyl groups on backbone
        backbone_oh = sum(1 for idx in path 
                         for n in mol.GetAtomWithIdx(idx).GetNeighbors()
                         if n.GetAtomicNum() == 8 and n.GetTotalNumHs() > 0)
        if backbone_oh < 1:
            return False, "Missing hydroxyl group(s) on glycerol backbone"

    except Exception as e:
        return False, f"Structural analysis failed: {str(e)}"

    # Optional: Check acyl chain length (at least 8 carbons)
    try:
        chain_start = mol.GetAtomWithIdx(ester_c).GetNeighbors()[0].GetIdx()
        chain = Chem.MolFromSmarts(f"[#{chain_start}](-[OX2])-*")
        chain_length = len(chain.GetSubstructMatch(mol)) - 1
        if chain_length < 8:
            return False, f"Acyl chain too short ({chain_length} carbons)"
    except:
        pass  # Skip if chain analysis fails

    return True, "Monoacylglycerol phosphate with proper 3-carbon backbone structure"