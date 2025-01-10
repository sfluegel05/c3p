"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester has an amine and ester group attached to the same carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern 1: Basic alpha-amino acid ester pattern
    # Excludes peptide bonds with negative lookahead
    pattern1 = Chem.MolFromSmarts("[NX3;!$(NC=O);!$(N=*);!$(N#*)][CH1,CH2][C](=O)O[CH2,CH3,CH]")
    
    # Pattern 2: N-substituted amino acid ester (including protected amines)
    pattern2 = Chem.MolFromSmarts("[#7;!$(NC=O);!$(N=*);!$(N#*)][CH1,CH2][C](=O)O[CH2,CH3,CH]")
    
    # Pattern 3: Cyclic amino acid esters (like proline derivatives)
    pattern3 = Chem.MolFromSmarts("[#7;R;!$(NC=O)][CH1,CH2;R][C](=O)O[CH2,CH3,CH]")
    
    # Pattern 4: Modified amino groups (azo compounds, etc.)
    pattern4 = Chem.MolFromSmarts("[#7;!$(NC=O)][CH1,CH2][C](=O)O[CH2,CH3,CH]")
    
    # Pattern 5: Special case for AMP esters
    pattern5 = Chem.MolFromSmarts("[NX3;!$(NC=O)][CH1][C](=O)O[C@@H]1[C@H](O)[C@H](O)[C@H](COP(=O)(O)O)O1")

    patterns = [
        (pattern1, "primary alpha-amino acid ester"),
        (pattern2, "N-substituted amino acid ester"),
        (pattern3, "cyclic amino acid ester"),
        (pattern4, "modified amino acid ester"),
        (pattern5, "nucleotide amino acid ester")
    ]
    
    for pattern, pattern_name in patterns:
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                n_atom = mol.GetAtomWithIdx(match[0])
                c_atom = mol.GetAtomWithIdx(match[1])
                
                # Verify nitrogen properties
                if n_atom.GetAtomicNum() != 7:
                    continue
                
                # Check for valid ester group
                ester_env = Chem.MolFromSmarts("[C](=O)O[CH2,CH3,CH]")
                if not mol.HasSubstructMatch(ester_env):
                    continue
                
                # Exclude cases where the ester is part of a macrocyclic peptide
                peptide_ring = Chem.MolFromSmarts("[C](=O)N[C]([C])C(=O)O[C]1[C][C]1")
                if mol.HasSubstructMatch(peptide_ring):
                    continue
                
                # Check alpha carbon environment
                alpha_c_neighbors = [atom.GetAtomicNum() for atom in c_atom.GetNeighbors()]
                if 7 in alpha_c_neighbors and 6 in alpha_c_neighbors:
                    # Additional check for azo compounds
                    if n_atom.GetIsAromatic() or "=" in Chem.MolToSmiles(mol, kekuleSmiles=True):
                        if pattern_name == "modified amino acid ester":
                            return True, f"Contains {pattern_name} pattern with modified amino group"
                    return True, f"Contains {pattern_name} pattern"
    
    return False, "No valid alpha-amino acid ester pattern found"