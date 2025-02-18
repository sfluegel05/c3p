"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:62705 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol (benzene ring with hydroxyl group) with exactly one nitro substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all benzene rings (6-membered aromatic carbocycles)
    ring_info = mol.GetRingInfo()
    benzene_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        all_carbon_aromatic = True
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                all_carbon_aromatic = False
                break
        if all_carbon_aromatic:
            benzene_rings.append(ring)
    
    # Precompute nitro group parent atoms (atoms attached to nitro groups)
    nitro_parents = []
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')
    for match in mol.GetSubstructMatches(nitro_pattern):
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        parent = None
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() not in match:
                parent = neighbor.GetIdx()
                break
        if parent is not None:
            nitro_parents.append(parent)
    
    # Check each benzene ring for hydroxyl and nitro count
    for ring in benzene_rings:
        # Check for hydroxyl group (-OH) on the ring
        has_hydroxyl = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetTotalNumHs() >= 1:
                        has_hydroxyl = True
                        break
            if has_hydroxyl:
                break
        if not has_hydroxyl:
            continue
        
        # Count nitro groups attached to this ring
        nitro_count = sum(1 for parent in nitro_parents if parent in ring)
        if nitro_count == 1:
            return True, "Contains a benzene ring with a hydroxyl group and exactly one nitro substituent"
    
    return False, "No benzene ring with hydroxyl and exactly one nitro group found"