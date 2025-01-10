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

    # SMARTS patterns for cyclohexenone
    # Pattern looks for:
    # - 6-membered ring
    # - One ketone (C=O)
    # - One C=C double bond
    # - All other bonds single
    # Multiple patterns to catch different arrangements
    patterns = [
        # Pattern 1: C=C-C-C-C-C(=O)
        "[#6]1=[#6]-[#6]-[#6]-[#6]-[#6](=O)1",
        # Pattern 2: C-C=C-C-C-C(=O)
        "[#6]1-[#6]=[#6]-[#6]-[#6]-[#6](=O)1",
        # Pattern 3: C-C-C=C-C-C(=O)
        "[#6]1-[#6]-[#6]=[#6]-[#6]-[#6](=O)1",
    ]

    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is None:
            continue
            
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            ring_atoms = set(match)
            
            # Check if any atom in the match is aromatic
            is_aromatic = False
            for atom_idx in ring_atoms:
                if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                    is_aromatic = True
                    break
                    
            if is_aromatic:
                continue
                
            # Verify ring atoms are all carbon (except oxygen of ketone)
            all_carbon = True
            for atom_idx in ring_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8:  # Allow C and O only
                    all_carbon = False
                    break
                    
            if not all_carbon:
                continue
                
            # Count double bonds in ring to verify pattern
            double_bonds = 0
            ketone_bonds = 0
            for bond in mol.GetBonds():
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                    
                start_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                
                # Check if bond involves ring atoms
                if start_idx in ring_atoms or end_idx in ring_atoms:
                    # Check if it's a ketone bond
                    if (mol.GetAtomWithIdx(start_idx).GetAtomicNum() == 8 or 
                        mol.GetAtomWithIdx(end_idx).GetAtomicNum() == 8):
                        ketone_bonds += 1
                    else:
                        double_bonds += 1
                        
            # Verify we have exactly one C=C double bond and one C=O bond
            if double_bonds == 1 and ketone_bonds == 1:
                return True, "Contains cyclohexenone ring (6-membered ring with one ketone and one double bond)"
                
    return False, "No valid cyclohexenone pattern found"