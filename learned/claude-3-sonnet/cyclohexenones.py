"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: cyclohexenones
Definition: Any six-membered alicyclic ketone having one double bond in the ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule contains a cyclohexenone structure.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_cyclohexenone, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Kekulize the molecule to ensure proper double bond representation
    Chem.Kekulize(mol)
    
    # SMARTS pattern for cyclohexenone
    # [C^2] represents sp2 carbon
    # The pattern looks for:
    # - 6-membered ring
    # - One ketone (C=O)
    # - One additional double bond in ring
    # - All carbons must be aliphatic (no aromatic bonds)
    cyclohexenone_pattern = Chem.MolFromSmarts("[C^2]1(=O)[C^2]=[C^2][C^2][C^2][C^2]1")
    
    # Find all matches
    matches = mol.GetSubstructMatches(cyclohexenone_pattern)
    
    if not matches:
        return False, "No cyclohexenone pattern found"
    
    # For each potential match, verify it's a valid cyclohexenone
    for match in matches:
        ring_atoms = set(match)
        
        # Check if any of the ring atoms are aromatic
        is_aromatic = False
        for atom_idx in ring_atoms:
            if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                is_aromatic = True
                break
                
        if is_aromatic:
            continue
            
        # Count double bonds in the ring
        double_bond_count = 0
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bond_count += 1
                    
        # Should have exactly 2 double bonds (one C=O and one C=C)
        if double_bond_count == 2:
            # Verify one is a ketone
            for atom_idx in ring_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 8:  # Oxygen
                    neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
                    if len(neighbors) == 1 and neighbors[0] == 6:  # Connected to exactly one carbon
                        return True, "Contains cyclohexenone ring (6-membered ring with one ketone and one double bond)"
                        
    return False, "No valid cyclohexenone pattern found (requires 6-membered ring with one ketone and one double bond)"