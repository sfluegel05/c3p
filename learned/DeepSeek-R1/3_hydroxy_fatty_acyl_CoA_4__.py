"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    Must have: 3-hydroxy acyl thioester, three phosphate groups with four deprotonated oxygens,
    and a linear fatty acid chain of at least 4 carbons.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for 3-hydroxy acyl thioester pattern (S-C(=O)-CH2-CH(OH)-*)
    thioester_pattern = Chem.MolFromSmarts("[S]-C(=O)-[CH2]-[CH](-O)")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing 3-hydroxy acyl thioester group"

    # Check acyl chain linearity and length (minimum 4 carbons total)
    min_chain_length = 4
    valid_chain_found = False
    
    for match in thioester_matches:
        try:
            s_idx, co_idx, c1_idx, c2_idx = match  # S, C=O, CH2, CH(OH)
            current_idx = c2_idx
            prev_idx = c1_idx
            chain_length = 4  # S-C(=O)-CH2-CH(OH) counts as 4 carbons

            while True:
                next_carbons = []
                for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
                    if (neighbor.GetAtomicNum() == 6 and 
                        neighbor.GetIdx() != prev_idx and 
                        neighbor.GetIdx() not in {s_idx, co_idx, c1_idx}):
                        next_carbons.append(neighbor.GetIdx())
                
                if len(next_carbons) > 1:
                    break  # Branch detected
                if not next_carbons:
                    break  # End of chain
                
                prev_idx = current_idx
                current_idx = next_carbons[0]
                chain_length += 1

            if chain_length >= min_chain_length and len(next_carbons) <= 1:
                valid_chain_found = True
                break
        except:
            continue

    if not valid_chain_found:
        return False, f"Acyl chain too short ({chain_length}) or branched"

    # Verify three phosphorus atoms (three phosphate groups)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 3:
        return False, f"Found {p_count} phosphorus atoms, need 3"

    # Count deprotonated oxygens (-O) attached to phosphorus
    o_minus_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                    o_minus_count += 1

    if o_minus_count != 4:
        return False, f"Found {o_minus_count} deprotonated phosphate oxygens, need 4"

    return True, "Contains 3-hydroxy acyl thioester, linear chain, 3 phosphates with 4 deprotonated oxygens"