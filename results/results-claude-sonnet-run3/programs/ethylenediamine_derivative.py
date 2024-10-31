from rdkit import Chem
from rdkit.Chem import AllChem

def is_ethylenediamine_derivative(smiles: str):
    """
    Determines if a molecule is an ethylenediamine derivative.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an ethylenediamine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Look for N-C-C-N pattern (ethylenediamine backbone)
    pattern = Chem.MolFromSmarts('[N]-[C]-[C]-[N]')
    if not mol.HasSubstructMatch(pattern):
        return False, "No ethylenediamine backbone (N-C-C-N) found"
        
    # Get matches of the pattern
    matches = mol.GetSubstructMatches(pattern)
    
    # For each match, verify:
    # 1. The carbons are sp3 hybridized (single bonds only)
    # 2. The nitrogens can have various substituents
    for match in matches:
        n1, c1, c2, n2 = match
        
        # Check carbon hybridization
        if mol.GetAtomWithIdx(c1).GetHybridization() != Chem.HybridizationType.SP3 or \
           mol.GetAtomWithIdx(c2).GetHybridization() != Chem.HybridizationType.SP3:
            continue
            
        # Check that carbons are connected by single bond
        bond = mol.GetBondBetweenAtoms(c1, c2)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
            
        # Verify both atoms are nitrogen
        if mol.GetAtomWithIdx(n1).GetSymbol() != 'N' or \
           mol.GetAtomWithIdx(n2).GetSymbol() != 'N':
            continue
            
        # Get substituents on nitrogens
        n1_atom = mol.GetAtomWithIdx(n1)
        n2_atom = mol.GetAtomWithIdx(n2)
        
        substituents = []
        for atom in [n1_atom, n2_atom]:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in [c1, c2]:
                    substituents.append(neighbor.GetSymbol())
                    
        if len(substituents) > 0:
            return True, f"Ethylenediamine derivative with substituents: {', '.join(set(substituents))}"
        else:
            return True, "Unsubstituted ethylenediamine"
            
    return False, "No valid ethylenediamine backbone found with correct bonding pattern"
# Pr=1.0
# Recall=1.0