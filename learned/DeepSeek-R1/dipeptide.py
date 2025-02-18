"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide contains two amino-acid residues connected by a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all amide bonds (C(=O)-N with single bond to N)
    # Adjusted pattern to include any N with three bonds (including proline-like cases)
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)-[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    valid_peptide_bonds = 0
    
    for match in amide_matches:
        c_atom_idx = match[0]
        n_atom_idx = match[1]
        
        # Get alpha carbon adjacent to carbonyl C (not part of the amide N)
        alpha1 = None
        c_atom = mol.GetAtomWithIdx(c_atom_idx)
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != n_atom_idx:
                alpha1 = neighbor
                break
        if not alpha1:
            continue  # Not a peptide bond
        
        # Get alpha carbon adjacent to amide N (not part of the carbonyl C)
        alpha2 = None
        n_atom = mol.GetAtomWithIdx(n_atom_idx)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != c_atom_idx:
                alpha2 = neighbor
                break
        if not alpha2:
            continue  # Not a peptide bond
        
        # Check substituents on alpha carbons (must have at least one side chain)
        # Check alpha1 (from C side) has at least one other neighbor besides carbonyl C
        alpha1_substituents = [n for n in alpha1.GetNeighbors() if n.GetIdx() != c_atom_idx]
        if not alpha1_substituents:
            continue
        
        # Check alpha2 (from N side) has at least one other neighbor besides amide N
        alpha2_substituents = [n for n in alpha2.GetNeighbors() if n.GetIdx() != n_atom_idx]
        if not alpha2_substituents:
            continue
        
        valid_peptide_bonds += 1
    
    if valid_peptide_bonds != 1:
        return False, f"Found {valid_peptide_bonds} valid peptide bonds, need exactly 1"
    
    return True, "Two amino acid residues connected by a peptide bond"