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
    
    # SMARTS pattern for 6-membered ring with a ketone
    # The * allows for any atom type in the ring (typically C but could be O,N etc)
    # The ketone oxygen must be double bonded to a ring atom
    pattern = Chem.MolFromSmarts("*1~*~*~*~*~*(=O)1")
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No 6-membered ring with ketone found"
        
    for match in matches:
        ring_atoms = set(match[:-1])  # Exclude the oxygen atom
        
        # Check if ring is aromatic
        is_aromatic = False
        for atom_idx in ring_atoms:
            if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                is_aromatic = True
                break
        if is_aromatic:
            continue
            
        # Count double bonds in ring (excluding the ketone)
        double_bond_count = 0
        ketone_carbon = None
        for bond in mol.GetBonds():
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            
            # Skip bonds not in ring
            if start_idx not in ring_atoms or end_idx not in ring_atoms:
                continue
                
            # Count double bonds
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bond_count += 1
                
        # Should have exactly 2 double bonds total (one ketone + one C=C)
        if double_bond_count == 2:
            return True, "Contains cyclohexenone ring (6-membered ring with one ketone and one double bond)"
            
    return False, "No valid cyclohexenone pattern found"