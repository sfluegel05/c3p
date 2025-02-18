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

    # Find all amide bonds (C(=O)-N)
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)-[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    valid_peptide_bonds = 0
    
    for match in amide_matches:
        c_atom_idx = match[0]
        n_atom_idx = match[1]
        
        # Get alpha carbon adjacent to carbonyl C (from first amino acid)
        alpha1 = None
        c_atom = mol.GetAtomWithIdx(c_atom_idx)
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != n_atom_idx:
                alpha1 = neighbor
                break
        if not alpha1:
            continue  # No alpha carbon on C side
        
        # Get alpha carbon adjacent to amide N (from second amino acid)
        alpha2 = None
        n_atom = mol.GetAtomWithIdx(n_atom_idx)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != c_atom_idx:
                alpha2 = neighbor
                break
        if not alpha2:
            continue  # No alpha carbon on N side
        
        valid_peptide_bonds += 1
    
    # Dipeptides must have exactly one peptide bond connecting two amino acids
    if valid_peptide_bonds != 1:
        return False, f"Found {valid_peptide_bonds} valid peptide bonds, need exactly 1"
    
    # Verify there are exactly two amino groups (one from each residue)
    # Check for at least two amine groups (alpha amines or side chains)
    amine_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check if it's an amine (NH2, NH, or part of a peptide bond)
            # Exclude the amide N in the peptide bond
            if atom.GetIdx() != n_atom_idx:
                # Check if it's part of an amino group (NH2 or NH)
                # Count primary, secondary amines (exclude amides)
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.SINGLE and atom.GetDegree() <= 3:
                        amine_count += 1
                        break
    
    if amine_count < 2:
        return False, f"Only {amine_count} amine groups found, need at least 2"
    
    return True, "Two amino acid residues connected by a peptide bond"