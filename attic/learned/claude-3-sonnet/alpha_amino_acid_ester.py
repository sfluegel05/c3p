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
    # Matches both primary and secondary amines
    pattern1 = Chem.MolFromSmarts("[NX3;!$(N=*);!$(N#*)][CH1,CH2][C](=O)O[C,H]")
    
    # Pattern 2: N-substituted amino acid ester pattern (including tertiary amines)
    pattern2 = Chem.MolFromSmarts("[#7;!$(N=*);!$(N#*)][CH1,CH2][C](=O)O[C,H]")
    
    # Pattern 3: Cyclic amino acid ester pattern (like proline derivatives)
    pattern3 = Chem.MolFromSmarts("[#7;R][CH1,CH2;R][C](=O)O[C,H]")
    
    # Pattern 4: Complex ring systems with embedded amino acid ester
    pattern4 = Chem.MolFromSmarts("[#7;R][CH1,CH2][C](=O)O[C,H]")
    
    patterns = [
        (pattern1, "basic alpha-amino acid ester"),
        (pattern2, "N-substituted amino acid ester"),
        (pattern3, "cyclic amino acid ester"),
        (pattern4, "ring-containing amino acid ester")
    ]
    
    for pattern, pattern_name in patterns:
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                n_atom = mol.GetAtomWithIdx(match[0])
                c_atom = mol.GetAtomWithIdx(match[1])
                
                # Verify nitrogen properties
                if n_atom.GetAtomicNum() != 7:  # Must be nitrogen
                    continue
                    
                # Check if carbon is really alpha to both N and C(=O)
                has_n = False
                has_ester = False
                for neighbor in c_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 7:  # Nitrogen
                        has_n = True
                    elif neighbor.GetAtomicNum() == 6:  # Carbon
                        for bond in neighbor.GetBonds():
                            other_atom = bond.GetOtherAtom(neighbor)
                            if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                                has_ester = True
                
                if has_n and has_ester:
                    # Additional check for valid ester group
                    ester_env = Chem.MolFromSmarts("[C](=O)O[C,H]")
                    if mol.HasSubstructMatch(ester_env):
                        return True, f"Contains {pattern_name} pattern"
    
    return False, "No valid alpha-amino acid ester pattern found"