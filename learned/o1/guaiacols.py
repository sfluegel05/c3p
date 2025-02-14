"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    
    A guaiacol is defined as any phenol carrying an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify phenolic OH groups
    phenolic_OH = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            if atom.GetDegree() == 1:  # Connected to one atom
                if atom.GetTotalNumHs(includeNeighbors=True) == 1:  # Has one hydrogen (OH group)
                    neighbor = atom.GetNeighbors()[0]
                    if neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6:
                        phenolic_OH.append((atom.GetIdx(), neighbor.GetIdx()))
    
    if not phenolic_OH:
        return False, "No phenolic hydroxyl group found"
    
    # For each phenolic OH group, check for ortho methoxy group
    for oh_idx, phenol_c_idx in phenolic_OH:
        phenol_c = mol.GetAtomWithIdx(phenol_c_idx)
        # Get neighboring carbons in the aromatic ring (ortho positions)
        ortho_carbons = []
        for neighbor in phenol_c.GetNeighbors():
            if neighbor.GetIdx() != oh_idx and neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6:
                ortho_carbons.append(neighbor)
        
        # Check if any ortho carbon has a methoxy group
        for ortho_c in ortho_carbons:
            for nbr in ortho_c.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 2:
                    # Oxygen atom connected to two atoms (methoxy oxygen)
                    attached_carbons = [atom.GetAtomicNum() for atom in nbr.GetNeighbors() if atom.GetAtomicNum() == 6]
                    attached_hydrogens = [atom.GetAtomicNum() for atom in nbr.GetNeighbors() if atom.GetAtomicNum() == 1]
                    if len(attached_carbons) == 1 and len(attached_hydrogens) == 0:
                        # Oxygen connected to one carbon and no hydrogens
                        methoxy_c = [atom for atom in nbr.GetNeighbors() if atom.GetAtomicNum() == 6 and atom.GetIdx() != ortho_c.GetIdx()][0]
                        if methoxy_c.GetDegree() == 1 and methoxy_c.GetTotalNumHs(includeNeighbors=True) == 3:
                            # Carbon is a methyl group (CH3)
                            return True, "Molecule is a guaiacol (phenolic OH with ortho methoxy group)"
    
    return False, "No ortho methoxy group found adjacent to phenolic hydroxyl group"