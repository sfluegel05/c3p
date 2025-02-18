"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidines have the structure RC(=NR)NR2 (a carbon double-bonded to an NR group and single-bonded to an NR2 group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern:
    # [CX3;!a] = Non-aromatic carbon with three bonds (double + two single)
    # =[NX2]   = Double bond to nitrogen with two connections (C and R group)
    # [NX3]    = Single bond to nitrogen with three connections (C and two R groups)
    pattern = Chem.MolFromSmarts("[CX3;!a](=[NX2])[NX3]")
    
    # Check for matches and verify nitrogen environments
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        # Additional validation to exclude cases where nitrogens are part of aromatic systems
        for match in matches:
            c_idx, n1_idx, n2_idx = match  # Assuming pattern order is C-N=N
            c_atom = mol.GetAtomWithIdx(c_idx)
            n1_atom = mol.GetAtomWithIdx(n1_idx)
            n2_atom = mol.GetAtomWithIdx(n2_idx)
            
            # Ensure none of the atoms are in aromatic rings
            if (c_atom.GetIsAromatic() or 
                n1_atom.GetIsAromatic() or 
                n2_atom.GetIsAromatic()):
                continue
                
            # Check hybridization states
            if (c_atom.GetHybridization() != Chem.HybridizationType.SP2 or
                n1_atom.GetHybridization() != Chem.HybridizationType.SP2 or
                n2_atom.GetHybridization() != Chem.HybridizationType.SP3):
                continue
                
            return True, "Contains carboxamidine group (RC(=NR)NR2)"
    
    return False, "No carboxamidine group detected"