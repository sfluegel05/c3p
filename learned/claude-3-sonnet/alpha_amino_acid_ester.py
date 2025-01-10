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

    # Core patterns for alpha-amino acid esters
    patterns = [
        # Primary alpha-amino acid ester
        (Chem.MolFromSmarts("[NX3;!$(NC=O);!$(N=*);!$(N#*)][CX4;!$(C=*);!$(C#*)][CX3](=O)[OX2][C,H1]"), 
         "primary alpha-amino acid ester"),
        
        # N-substituted alpha-amino acid ester (including protected amines)
        (Chem.MolFromSmarts("[#7;!$(NC=O);!$(N=*);!$(N#*);!$(N@[a])]([#6,#1])[CX4;!$(C=*);!$(C#*)][CX3](=O)[OX2][C,H1]"),
         "N-substituted alpha-amino acid ester"),
        
        # Cyclic amino acid ester (proline-type)
        (Chem.MolFromSmarts("[#7;R;!$(NC=O);!$(N=*);!$(N#*)][CX4;R;!$(C=*);!$(C#*)][CX3](=O)[OX2][C,H1]"),
         "cyclic alpha-amino acid ester"),
        
        # Special case for nucleotide amino acid esters
        (Chem.MolFromSmarts("[NX3;!$(NC=O)][CX4;!$(C=*);!$(C#*)][CX3](=O)[OX2][C@@H]1[C@H](O)[C@H](O)[C@H](COP(=O)(O)O)O1"),
         "nucleotide amino acid ester")
    ]

    for pattern, pattern_name in patterns:
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                # Get key atoms
                n_atom = mol.GetAtomWithIdx(match[0])
                c_alpha = mol.GetAtomWithIdx(match[1])
                c_carbonyl = mol.GetAtomWithIdx(match[2])
                
                # Verify nitrogen properties
                if n_atom.GetAtomicNum() != 7:
                    continue
                    
                # Check that nitrogen is not part of aromatic system
                if n_atom.GetIsAromatic():
                    continue
                    
                # Verify alpha carbon properties
                if c_alpha.GetAtomicNum() != 6:
                    continue
                    
                # Check that alpha carbon is sp3
                if c_alpha.GetHybridization() != Chem.HybridizationType.SP3:
                    continue
                    
                # Verify carbonyl carbon properties
                if c_carbonyl.GetAtomicNum() != 6:
                    continue
                    
                # Additional checks for cyclic systems
                if pattern_name == "cyclic alpha-amino acid ester":
                    ring_info = mol.GetRingInfo()
                    if not ring_info.IsAtomInRingOfSize(match[0], 5) and not ring_info.IsAtomInRingOfSize(match[0], 6):
                        continue
                
                # Exclude cases where the amine is part of complex ring systems
                complex_ring = Chem.MolFromSmarts("[#7]1@[#6]@[#6]@[#6]@[#6]@[#6]@[#6]@1")
                if mol.HasSubstructMatch(complex_ring):
                    if not any(match[0] in m for m in mol.GetSubstructMatches(complex_ring)):
                        return True, f"Contains {pattern_name} pattern"
                else:
                    return True, f"Contains {pattern_name} pattern"
                
    return False, "No valid alpha-amino acid ester pattern found"