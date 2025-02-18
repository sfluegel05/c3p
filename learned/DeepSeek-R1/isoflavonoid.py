"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran (chromen-4-one) with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the chromen-4-one core with labeled C3 atom (position 3)
    core_smarts = Chem.MolFromSmarts('c1ccc2c(c1)oc(=O)c([*:3])cc2')
    if not core_smarts:
        return False, "Invalid core SMARTS pattern"
    
    matches = mol.GetSubstructMatches(core_smarts)
    if not matches:
        return False, "No chromen-4-one core found"
    
    # Locate the C3 atom in the core pattern
    c3_core_idx = None
    for atom in core_smarts.GetAtoms():
        if atom.HasProp('molAtomMapNumber') and atom.GetProp('molAtomMapNumber') == '3':
            c3_core_idx = atom.GetIdx()
            break
    if c3_core_idx is None:
        return False, "C3 position not identified in core"
    
    # Check each core instance for an aryl substituent at C3
    for match in matches:
        c3_mol_idx = match[c3_core_idx]
        c3_atom = mol.GetAtomWithIdx(c3_mol_idx)
        
        for neighbor in c3_atom.GetNeighbors():
            if neighbor.GetIdx() in match:  # Skip atoms in the core
                continue
            
            # Check if neighbor is part of an aromatic ring (>=6 atoms)
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if neighbor.GetIdx() in ring:
                    if len(ring) >= 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                        return True, "Chromen-4-one core with aryl group at position 3"
    
    return False, "No aryl substituent at position 3"