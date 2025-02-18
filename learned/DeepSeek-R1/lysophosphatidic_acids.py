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
    - Hydroxyl group on the remaining backbone carbon
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Find exactly one ester group (O-C=O)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]C(=O)"))
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups (need 1)"
    
    # Get glycerol carbon attached to ester oxygen
    ester_o = ester_matches[0][0]
    ester_carbonyl = ester_matches[0][1]
    glycerol_c_ester = [n.GetIdx() for n in mol.GetAtomWithIdx(ester_o).GetNeighbors() 
                       if n.GetIdx() != ester_carbonyl and n.GetAtomicNum() == 6][0]

    # Find phosphate group (P=O connected to carbon via oxygen)
    phosphate_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[PX4](=O)(O)(O)OC"))
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Get glycerol carbon attached to phosphate
    phosphate_p = phosphate_matches[0][0]
    phosphate_o = [n for n in mol.GetAtomWithIdx(phosphate_p).GetNeighbors() 
                  if n.GetAtomicNum() == 8 and 
                  any(nn.GetAtomicNum() == 6 for nn in n.GetNeighbors())][0]
    glycerol_c_phosphate = [n.GetIdx() for n in phosphate_o.GetNeighbors() 
                           if n.GetAtomicNum() == 6][0]

    # Check backbone connectivity
    path = Chem.GetShortestPath(mol, glycerol_c_ester, glycerol_c_phosphate)
    if len(path) not in [2, 3]:  # Allow adjacent or opposite carbons
        return False, f"Invalid backbone path length ({len(path)}) between ester/phosphate"
    
    # Verify 3-carbon backbone
    backbone = set(path)
    for atom_idx in path:
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6:
            return False, "Non-carbon in backbone"
    
    # Find third carbon with hydroxyl
    hydroxyl_found = False
    for carbon in backbone:
        for neighbor in mol.GetAtomWithIdx(carbon).GetNeighbors():
            if (neighbor.GetAtomicNum() == 8 and 
                neighbor.GetTotalNumHs() > 0 and 
                mol.GetBondBetweenAtoms(carbon, neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE):
                hydroxyl_found = True
    if not hydroxyl_found:
        return False, "Missing hydroxyl on backbone"

    # Optional acyl chain length check
    try:
        chain_start = ester_carbonyl
        chain_length = 1  # Starts at carbonyl carbon
        while True:
            next_c = [n for n in mol.GetAtomWithIdx(chain_start).GetNeighbors() 
                     if n.GetAtomicNum() == 6 and n.GetIdx() != ester_o]
            if not next_c:
                break
            chain_start = next_c[0].GetIdx()
            chain_length += 1
        if chain_length < 8:
            return False, f"Acyl chain too short ({chain_length} carbons)"
    except:
        pass

    return True, "Monoacylglycerol phosphate with valid 3-carbon backbone structure"