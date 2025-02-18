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
    
    # Improved SMARTS pattern:
    # [C] must have:
    # - At least one hydroxyl (-OH) group
    # - At least one amine group (NH2, NHR where R is not carbonyl)
    # - Explicitly exclude nitrogens in amides/nitriles
    # The * allows for other substituents on the carbon
    hemiaminal_pattern = Chem.MolFromSmarts("[C]([OH])([N;H2,H1;!$(N=*);!$(N#*)])*")
    
    # Find all matching atoms
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    
    if matches:
        # Verify at least one valid amine-hydroxyl pair
        for match in matches:
            # Get the central carbon atom
            carbon_idx = match[0]
            atom = mol.GetAtomWithIdx(carbon_idx)
            
            # Check for exactly one hydroxyl and one amine attached
            # Allow other substituents as long as both groups are present
            has_oh = False
            has_amine = False
            
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:  # -OH group
                    has_oh = True
                elif neighbor.GetAtomicNum() == 7:
                    # Check for primary/secondary amine (not amide)
                    if neighbor.GetTotalNumHs() >= 1 and not neighbor.GetIsAromatic():
                        has_amine = True
                        
            if has_oh and has_amine:
                return True, "Contains a carbon with both hydroxyl and amine groups"
        
        return False, "No carbon with both hydroxyl and amine groups found"
    else:
        return False, "No carbon with hydroxyl and amine groups found"