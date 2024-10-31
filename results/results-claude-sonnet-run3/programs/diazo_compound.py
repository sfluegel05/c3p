from rdkit import Chem
from rdkit.Chem import AllChem

def is_diazo_compound(smiles: str):
    """
    Determines if a molecule is a diazo compound (contains =N2 group attached to single carbon).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a diazo compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for diazo group pattern using SMARTS
    # [C]=[N+]=[N-] or [C]=[N]=[N] or [C]=N#N patterns
    diazo_patterns = [
        '[C]=[N+]=[N-]',
        '[C]=[N]=[N]',
        '[C]=N#N'
    ]
    
    for pattern in diazo_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if matches:
            # Get the carbon atom that's part of the diazo group
            for match in matches:
                carbon_idx = match[0]
                carbon = mol.GetAtomWithIdx(carbon_idx)
                
                # Check that carbon is attached to exactly one diazo group
                n_diazo = 0
                for neighbor in carbon.GetNeighbors():
                    if neighbor.GetSymbol() == 'N':
                        n2_group = False
                        for n2_neighbor in neighbor.GetNeighbors():
                            if n2_neighbor.GetSymbol() == 'N':
                                n2_group = True
                                break
                        if n2_group:
                            n_diazo += 1
                
                if n_diazo == 1:
                    return True, f"Contains diazo group matching pattern {pattern}"
                    
    return False, "No diazo group found"
# Pr=1.0
# Recall=1.0