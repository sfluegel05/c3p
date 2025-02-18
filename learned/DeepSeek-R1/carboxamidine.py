"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for carboxamidine core: C(=N)-N with appropriate substituents
    pattern = Chem.MolFromSmarts("[C;!a](=[N;D2])[N;D3]")
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No carboxamidine group detected"
    
    for match in matches:
        c_idx, n1_idx, n2_idx = match  # Assuming pattern order is C, N1, N2
        
        c_atom = mol.GetAtomWithIdx(c_idx)
        n1_atom = mol.GetAtomWithIdx(n1_idx)
        n2_atom = mol.GetAtomWithIdx(n2_idx)
        
        # Check none of the atoms are aromatic
        if (c_atom.GetIsAromatic() or 
            n1_atom.GetIsAromatic() or 
            n2_atom.GetIsAromatic()):
            continue
        
        # Check degrees are correct (redundant due to SMARTS, but safe)
        if n1_atom.GetDegree() != 2 or n2_atom.GetDegree() != 3:
            continue
        
        return True, "Contains carboxamidine group (RC(=NR)NR2)"
    
    return False, "No valid carboxamidine group found"