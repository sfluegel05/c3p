"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:50047 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has an amino group and a hydroxy group attached to the same carbon atom.
    Validates both groups are on the same carbon and the nitrogen is in amine form (not amide/nitrile).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to find carbons with both hydroxyl and amine groups
    # [C] must have:
    # - At least one hydroxyl (-O-H)
    # - At least one amine (NH2, NHR where R is not part of amide/nitrile)
    # Excludes amides (N-C=O), aromatic nitrogens, and nitriles
    hemiaminal_pattern = Chem.MolFromSmarts("[C][OH]-&@[C](-[N;H2,H1;!$(N-C=[O,S,N]);!$(N#*)])")
    
    # Alternative approach: Check for any carbon with both OH and NH groups attached
    # This pattern matches a carbon connected to both OH and NH (primary/secondary amine)
    hemiaminal_pattern = Chem.MolFromSmarts("[C]([OH])([N;H2,H1;!$(N-C=[O,S,N]);!$(N#*)])")
    
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    
    if not matches:
        return False, "No carbon with hydroxyl and amine groups found"
    
    # Verify each matched carbon has at least one OH and one amine
    for match in matches:
        carbon_idx = match[0]
        atom = mol.GetAtomWithIdx(carbon_idx)
        has_oh = False
        has_amine = False
        
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                has_oh = True
            elif neighbor.GetAtomicNum() == 7:
                # Check if it's a non-amide, non-nitrile amine
                if neighbor.GetTotalNumHs() >= 1 and not any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor.GetBonds()):
                    has_amine = True
        
        if has_oh and has_amine:
            return True, "Contains a carbon with both hydroxyl and amine groups"
    
    return False, "No carbon with both hydroxyl and amine groups found"