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
    - Glycerol backbone (3 carbons in chain)
    - Exactly one ester-linked acyl chain
    - Phosphate group attached to backbone
    - Hydroxyl group on remaining backbone carbon
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Find exactly one ester group (O-C=O connected to carbon chain)
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups (need 1)"
    
    # Get glycerol carbon connected to ester oxygen
    ester_o = ester_matches[0][0]
    ester_carbon = [n for n in mol.GetAtomWithIdx(ester_o).GetNeighbors() 
                   if n.GetIdx() != ester_matches[0][1] and n.GetAtomicNum() == 6]
    if not ester_carbon:
        return False, "Ester oxygen not connected to backbone"
    glycerol_c_ester = ester_carbon[0].GetIdx()

    # Find phosphate group (P=O with two OH and connected to backbone)
    phosphate_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[PX4](=O)(O)(O)OC"))
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Get glycerol carbon connected to phosphate
    phosphate_p = phosphate_matches[0][0]
    phosphate_o = next((n for n in mol.GetAtomWithIdx(phosphate_p).GetNeighbors() 
                       if n.GetAtomicNum() == 8 and 
                       any(nn.GetAtomicNum() == 6 for nn in n.GetNeighbors())), None)
    if not phosphate_o:
        return False, "Phosphate not connected to backbone"
    
    glycerol_c_phosphate = [n.GetIdx() for n in phosphate_o.GetNeighbors() if n.GetAtomicNum() == 6]
    if not glycerol_c_phosphate:
        return False, "Phosphate oxygen not connected to carbon"
    glycerol_c_phosphate = glycerol_c_phosphate[0]

    # Verify 3-carbon backbone connectivity
    try:
        path = Chem.GetShortestPath(mol, glycerol_c_ester, glycerol_c_phosphate)
    except:
        return False, "No path between ester and phosphate carbons"
    
    backbone = set(path)
    # Check for third carbon in glycerol chain
    if len(backbone) == 2:
        # Look for adjacent carbon to complete triplet
        adj1 = [n.GetIdx() for n in mol.GetAtomWithIdx(glycerol_c_ester).GetNeighbors() 
               if n.GetAtomicNum() == 6 and n.GetIdx() not in backbone]
        adj2 = [n.GetIdx() for n in mol.GetAtomWithIdx(glycerol_c_phosphate).GetNeighbors() 
               if n.GetAtomicNum() == 6 and n.GetIdx() not in backbone]
        if adj1:
            backbone.add(adj1[0])
        elif adj2:
            backbone.add(adj2[0])
        else:
            return False, "Cannot find third backbone carbon"
    
    if len(backbone) != 3:
        return False, f"Backbone has {len(backbone)} carbons (need 3)"
    
    # Check all backbone atoms are carbons
    for c in backbone:
        if mol.GetAtomWithIdx(c).GetAtomicNum() != 6:
            return False, "Non-carbon in backbone"

    # Verify hydroxyl group on third carbon
    hydroxyl_found = False
    for c in backbone:
        if c in {glycerol_c_ester, glycerol_c_phosphate}:
            continue
        for bond in mol.GetAtomWithIdx(c).GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtomIdx(c)
                if (mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 8 and 
                    mol.GetAtomWithIdx(neighbor).GetTotalNumHs() > 0):
                    hydroxyl_found = True
    if not hydroxyl_found:
        return False, "Missing hydroxyl on third carbon"

    return True, "Monoacylglycerol phosphate with valid 3-carbon backbone"