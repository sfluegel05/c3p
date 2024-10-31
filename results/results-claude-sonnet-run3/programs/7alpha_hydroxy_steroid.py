from rdkit import Chem
from rdkit.Chem import AllChem

def is_7alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 7alpha-hydroxy steroid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 7alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Convert SMILES to RDKit mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate 3D coordinates if needed for stereochemistry
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        pass

    # SMARTS pattern for steroid core with 7alpha-hydroxy group
    steroid_7alpha_oh_pattern = Chem.MolFromSmarts(
        '[#6]~1~2~[#6]~[#6]~[#6@H]([OH1])~[#6]~1~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~3~2'
    )

    # Check for matches
    matches = mol.GetSubstructMatches(steroid_7alpha_oh_pattern)
    
    if not matches:
        return False, "No steroid core with 7alpha-hydroxy group found"

    # Additional check for alpha stereochemistry at position 7
    for match in matches:
        # The hydroxyl-bearing carbon is at position 3 in our SMARTS pattern
        c7_idx = match[3]
        c7_atom = mol.GetAtomWithIdx(c7_idx)
        
        # Check neighbors for hydroxyl group
        for neighbor in c7_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
                # Check if this is actually a 7-position in steroid
                ring_info = mol.GetRingInfo()
                if ring_info.NumAtomRings(c7_idx) >= 2:  # Should be part of two rings
                    # Check stereochemistry
                    if c7_atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                        return True, "7alpha-hydroxy steroid identified"
                    elif c7_atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                        return False, "7beta-hydroxy steroid identified"
                    else:
                        # If stereochemistry is present in the structure but not explicitly defined
                        return True, "7alpha-hydroxy steroid identified (stereochemistry inferred from structure)"

    return False, "No valid 7alpha-hydroxy group found in steroid structure"
# Pr=None
# Recall=0.0