"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:31011 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester (laurate ester) based on its SMILES string.
    A dodecanoate ester contains an ester group where the acyl chain is from lauric acid (12 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for ester group O-C(=O)-C (acyl chain start)
    ester_pattern = Chem.MolFromSmarts('[OX2]-C(=O)-[!O]')  # O connected to carbonyl, then non-O atom
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not ester_matches:
        return False, "No ester group found"
    
    # Check each ester group for C12 acyl chain
    for match in ester_matches:
        if len(match) != 3:
            continue  # Skip invalid matches (shouldn't occur with this pattern)
        o_idx, co_idx, acyl_start_idx = match
        
        acyl_start_atom = mol.GetAtomWithIdx(acyl_start_idx)
        if acyl_start_atom.GetAtomicNum() != 6:
            continue  # Acyl chain must start with carbon
        
        # Traverse the acyl chain to count length
        current = acyl_start_atom
        visited = {co_idx, current.GetIdx()}
        count = 1  # Starting at first carbon after carbonyl
        
        while True:
            next_carbons = []
            for neighbor in current.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx in visited:
                    continue
                if neighbor.GetAtomicNum() == 6 and neighbor.GetTotalNumHs() >= 2:  # Chain carbons
                    next_carbons.append(neighbor)
            
            # Check for unbranched chain
            if len(next_carbons) != 1:
                break  # Branch or end of chain
            current = next_carbons[0]
            visited.add(current.GetIdx())
            count += 1
            
            if count > 11:  # Safety break for long chains
                break
        
        if count == 11:  # 11 carbons after carbonyl (total 12 with carbonyl)
            return True, "Contains a dodecanoyl (C12) ester group"
    
    return False, "No dodecanoyl ester group found"