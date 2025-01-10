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
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH,O-]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Look for N-acetyl group (CH3-C(=O)-N)
    n_acetyl = Chem.MolFromSmarts('[CH3]C(=O)[NH0]')  # NH0 means nitrogen with no hydrogens
    if not mol.HasSubstructMatch(n_acetyl):
        return False, "No N-acetyl group found"

    # More general pattern for N-acetyl amino acid core
    # Matches: CH3-C(=O)-N-C-C(=O)O where the alpha carbon can have any substitution
    core_pattern = Chem.MolFromSmarts('[CH3]C(=O)[NH0]C([*,H])[*]C(=O)[OH,O-]')
    
    # Alternative pattern for cyclic amino acids (like proline derivatives)
    cyclic_pattern = Chem.MolFromSmarts('[CH3]C(=O)N1[CH2,CH]CC[CH2,CH]C1C(=O)[OH,O-]')
    
    if not (mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(cyclic_pattern)):
        return False, "No N-acetyl amino acid core structure found"

    # Count carbons and nitrogens to avoid false positives
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 4:  # Minimum carbons needed (CH3CO-NH-CH-COOH)
        return False, "Too few carbons for N-acetyl amino acid"
    
    # Verify that the N-acetyl group is connected to the amino acid part
    # by checking if there's a path between the N-acetyl nitrogen and the carboxylic acid
    n_acetyl_matches = mol.GetSubstructMatches(n_acetyl)
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid)
    
    for n_match in n_acetyl_matches:
        n_idx = n_match[2]  # Index of nitrogen atom in N-acetyl match
        for c_match in carboxyl_matches:
            c_idx = c_match[0]  # Index of carbon atom in carboxyl match
            # Check if there's a path of appropriate length between N and COOH
            path = Chem.GetShortestPath(mol, n_idx, c_idx)
            if path and 2 <= len(path) <= 6:  # Path length should be reasonable
                return True, "Contains N-acetyl group properly connected to amino acid core"
                
    return False, "N-acetyl group not properly connected to amino acid core"