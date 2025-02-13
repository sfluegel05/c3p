"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl-amino acid
Definition: An N-acyl-amino acid that has acetyl as the acyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (including deprotonated form)
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1,OX1-]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Look for N-acetyl group (CH3-C(=O)-N)
    n_acetyl = Chem.MolFromSmarts('[CH3][CX3](=O)[NX3]')
    if not mol.HasSubstructMatch(n_acetyl):
        return False, "No N-acetyl group found"

    # Look for amino acid core with N-acetyl
    # More flexible pattern that allows different substitutions on alpha carbon
    amino_acid_core = Chem.MolFromSmarts('[CH3][CX3](=O)[NX3][CX4][CX3](=O)[OX2H1,OX1-]')
    if not mol.HasSubstructMatch(amino_acid_core):
        return False, "No N-acetylated amino acid core structure found"

    # Additional check for cyclic amino acids (like proline derivatives)
    cyclic_amino_acid = Chem.MolFromSmarts('[CH3]C(=O)N1[CH2,CH]CC[CH2,CH][C@H]1C(=O)[OH]')
    
    # Count number of acetyl groups
    n_acetyl_matches = len(mol.GetSubstructMatches(n_acetyl))
    if n_acetyl_matches > 2:
        return False, f"Too many N-acetyl groups ({n_acetyl_matches})"

    # Verify connectivity: N-acetyl must be on amino acid nitrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            neighbors = atom.GetNeighbors()
            acetyl_c = False
            alpha_c = False
            for neighbor in neighbors:
                # Check if neighbor is carbonyl carbon of acetyl
                if neighbor.GetAtomicNum() == 6:
                    for n2 in neighbor.GetNeighbors():
                        if n2.GetAtomicNum() == 8 and n2.GetDegree() == 1:  # C=O
                            for n3 in neighbor.GetNeighbors():
                                if n3.GetAtomicNum() == 6 and n3.GetDegree() == 4:  # CH3
                                    acetyl_c = True
                # Check if neighbor is alpha carbon (connected to COOH)
                if neighbor.GetAtomicNum() == 6:
                    for n2 in neighbor.GetNeighbors():
                        if n2.GetAtomicNum() == 6:
                            for n3 in n2.GetNeighbors():
                                if n3.GetAtomicNum() == 8:
                                    alpha_c = True

            if acetyl_c and alpha_c:
                return True, "Contains N-acetyl group attached to amino acid core"
            
    # Check for cyclic amino acids separately
    if mol.HasSubstructMatch(cyclic_amino_acid):
        return True, "Contains N-acetyl group attached to cyclic amino acid core"

    return False, "Structure does not match N-acetyl-amino acid pattern"