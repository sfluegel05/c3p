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
    
    # Neutralize the molecule to handle different protonation states
    uncharger = Chem.AllChem.UnchargeModel()
    mol = uncharger.uncharge(mol)

    # Look for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1,OX1-]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Look for N-acetyl group (CH3-C(=O)-N)
    n_acetyl = Chem.MolFromSmarts('[CH3][CX3](=O)[NX3]')
    if not mol.HasSubstructMatch(n_acetyl):
        return False, "No N-acetyl group found"

    # Look for amino acid core with N-acetyl
    # [C] is alpha carbon, [N] is acetylated N, [C]=O is carboxylic acid
    amino_acid_core = Chem.MolFromSmarts('[CH3][CX3](=O)[NX3][CX4][CX3](=O)[OX2H1,OX1-]')
    if not mol.HasSubstructMatch(amino_acid_core):
        return False, "No N-acetylated amino acid core structure found"

    # Count number of N-acetyl groups
    n_acetyl_matches = len(mol.GetSubstructMatches(n_acetyl))
    if n_acetyl_matches > 2:
        return False, f"Too many N-acetyl groups ({n_acetyl_matches})"

    # Count carbons and oxygens to ensure reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbons for N-acetyl-amino acid"

    # Check that nitrogen is connected to both acetyl and alpha carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            neighbors = atom.GetNeighbors()
            acetyl_c = False
            alpha_c = False
            for neighbor in neighbors:
                # Check if neighbor is carbonyl carbon of acetyl
                if neighbor.GetAtomicNum() == 6 and any(n.GetAtomicNum() == 8 and n.GetDegree() == 1 
                                                      for n in neighbor.GetNeighbors()):
                    acetyl_c = True
                # Check if neighbor is alpha carbon (connected to COOH)
                if neighbor.GetAtomicNum() == 6 and any(n.GetAtomicNum() == 6 and 
                                                      any(nn.GetAtomicNum() == 8 for nn in n.GetNeighbors())
                                                      for n in neighbor.GetNeighbors()):
                    alpha_c = True
            if acetyl_c and alpha_c:
                return True, "Contains N-acetyl group attached to amino acid core"

    return False, "Structure does not match N-acetyl-amino acid pattern"