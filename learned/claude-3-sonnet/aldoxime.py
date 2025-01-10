"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:33566 aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    Aldoximes are compounds with the general structure H-CR=N-OH, where R is any carbon group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Main pattern for aldoxime group: CH=N-OH where:
    # - C has exactly one H
    # - C is not aromatic
    # - C is not part of a C=O group
    # - N is not charged
    # - O has exactly one H
    pattern = "[CH1X3!a;!$(C=O);!$(C[N+]);!$(C=[N+])]=[NX2v3;!+]-[OH1X2]"
    
    smarts = Chem.MolFromSmarts(pattern)
    if not mol.HasSubstructMatch(smarts):
        return False, "No aldoxime group (H-CR=N-OH) found"
        
    matches = mol.GetSubstructMatches(smarts)
    
    for match in matches:
        carbon_idx = match[0]
        nitrogen_idx = match[1]
        oxygen_idx = match[2]
        
        carbon = mol.GetAtomWithIdx(carbon_idx)
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        
        # Additional validations:
        
        # 1. Verify carbon has exactly one hydrogen
        if carbon.GetTotalNumHs() != 1:
            continue
            
        # 2. Verify carbon is not in a ring
        if carbon.IsInRing():
            ring_info = mol.GetRingInfo()
            if any(carbon_idx in ring for ring in ring_info.AtomRings()):
                continue
                
        # 3. Check that nitrogen has no other double bonds
        n_bonds = nitrogen.GetBonds()
        if sum(1 for b in n_bonds if b.GetBondType() == Chem.BondType.DOUBLE) > 1:
            continue
            
        # 4. Verify oxygen is only connected to nitrogen
        o_bonds = oxygen.GetBonds()
        if len(o_bonds) != 1 or oxygen.GetTotalNumHs() != 1:
            continue
            
        # 5. Exclude cases where carbon is part of activated systems
        carbon_neighbors = carbon.GetNeighbors()
        skip = False
        for neighbor in carbon_neighbors:
            if neighbor.GetAtomicNum() == 7 and neighbor != nitrogen:  # No other N attachments
                skip = True
                break
            if neighbor.GetIsAromatic():  # Not directly connected to aromatic systems
                skip = True
                break
        if skip:
            continue
            
        # If we pass all validations, it's an aldoxime
        return True, "Contains aldoxime group (H-CR=N-OH)"
        
    return False, "No aldoxime group (H-CR=N-OH) found"