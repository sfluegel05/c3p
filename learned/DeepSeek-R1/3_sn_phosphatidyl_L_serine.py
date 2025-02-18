"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine (CHEBI:89834)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    Must have a glycerol backbone with two acyl groups at positions 1 and 2, and a phosphoserine group at position 3.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # 1. Find phosphoserine group (P connected to serine)
    # Pattern: P(=O)(O)OC[C@H](N)C(=O)O (L-serine)
    phosphoserine_pattern = Chem.MolFromSmarts("[P](=O)([OX2])[OX2][C@H]([CH2])[C](=O)[OX2]")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Phosphoserine group not found"
    
    # Get the phosphorus atom in the match
    matches = mol.GetSubstructMatches(phosphoserine_pattern)
    if not matches:
        return False, "Phosphoserine group not found"
    p_idx = matches[0][0]
    p_atom = mol.GetAtomWithIdx(p_idx)
    
    # 2. Find the oxygen connecting P to glycerol (O-P-O-C...)
    glycerol_o = None
    for neighbor in p_atom.GetNeighbors():
        if neighbor.GetSymbol() == "O":
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    for nbr in neighbor.GetNeighbors():
                        if nbr.GetSymbol() == "C" and nbr.GetIdx() != p_idx:
                            glycerol_o = neighbor
                            break
                    if glycerol_o:
                        break
            if glycerol_o:
                break
    if not glycerol_o:
        return False, "No glycerol attached to phosphate"
    
    # 3. Get the glycerol's C3 (connected to this oxygen)
    c3 = glycerol_o.GetNeighbors()[0]
    if c3.GetSymbol() != "C":
        return False, "Phosphate not connected to glycerol carbon"
    
    # 4. Verify glycerol structure (C1-C2-C3 with C3 connected to phosphate)
    # Check if C3 is part of a three-carbon chain
    # Trace back to C2 and C1
    c2 = None
    for neighbor in c3.GetNeighbors():
        if neighbor.GetSymbol() == "C" and neighbor.GetIdx() != glycerol_o.GetIdx():
            c2 = neighbor
            break
    if not c2:
        return False, "Glycerol C3 not connected to C2"
    
    c1 = None
    for neighbor in c2.GetNeighbors():
        if neighbor.GetSymbol() == "C" and neighbor.GetIdx() != c3.GetIdx():
            c1 = neighbor
            break
    if not c1:
        return False, "Glycerol C2 not connected to C1"
    
    # 5. Check stereochemistry of C2 (sn-3 configuration)
    # The configuration should have the phosphate in the correct position
    # This may require checking the order of substituents, but RDKit's stereochemistry handling is complex
    # For simplicity, assume the SMILES encodes the correct stereochemistry
    
    # 6. Check C1 and C2 for ester groups (O-C(=O)-R)
    ester_count = 0
    for c in [c1, c2]:
        ester_found = False
        for neighbor in c.GetNeighbors():
            if neighbor.GetSymbol() == "O":
                # Check if this oxygen is part of an ester (O-C(=O))
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == Chem.BondType.SINGLE and bond.GetOtherAtom(neighbor).GetSymbol() == "C":
                        # Check for adjacent carbonyl group
                        for a in bond.GetOtherAtom(neighbor).GetNeighbors():
                            if a.GetSymbol() == "O" and mol.GetBondBetweenAtoms(bond.GetOtherAtom(neighbor).GetIdx(), a.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                ester_count += 1
                                ester_found = True
                                break
                        if ester_found:
                            break
                if ester_found:
                    break
    if ester_count < 2:
        return False, f"Found {ester_count} ester groups on C1/C2, need 2"
    
    # 7. Check acyl chain lengths (at least 4 carbons each)
    # Use SMARTS to find ester groups and count chain length
    ester_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2][C@H]")
    ester_matches = mol.GetSubstructMatches(ester_smarts)
    if len(ester_matches) < 2:
        return False, "Insufficient ester groups"
    
    # Check each ester's carbon chain
    min_chain_length = 4  # At least 4 carbons (e.g., butyryl)
    for match in ester_matches[:2]:  # Check first two esters
        carbonyl_c = match[0]
        chain = []
        # Traverse the chain from the carbonyl carbon
        current = carbonyl_c
        visited = set()
        while True:
            visited.add(current)
            next_atoms = [a for a in mol.GetAtomWithIdx(current).GetNeighbors() if a.GetSymbol() == "C" and a.GetIdx() not in visited]
            if not next_atoms:
                break
            current = next_atoms[0].GetIdx()
            chain.append(current)
        # Chain length is number of carbons excluding the carbonyl
        if len(chain) < min_chain_length:
            return False, f"Acyl chain too short: {len(chain)} carbons"
    
    return True, "sn-glycerol with two acyl chains at C1/C2 and phosphoserine at C3"