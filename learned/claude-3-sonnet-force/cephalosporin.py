"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: CHEBI:35469 cephalosporin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are a class of beta-lactam antibiotics with a 6-membered dihydrothiazine side ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for beta-lactam ring pattern (N-C-C-N-C=O)
    beta_lactam_pattern = Chem.MolFromSmarts("[N,O]=[C-]1[N,O]C(=O)C1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Look for dihydrothiazine ring pattern (6-membered ring with N, S, and C=C)
    dihydrothiazine_pattern = Chem.MolFromSmarts("C1NC(=O)C(=C)S1")
    if not mol.HasSubstructMatch(dihydrothiazine_pattern):
        return False, "No dihydrothiazine ring found"

    # Check for common structural features
    has_amino_acid_side_chain = mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]C(=O)[CH2]"))
    has_heteroaromatic_ring = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in ["c1ncnc1", "c1ccncc1", "c1cncnc1"])
    has_oxyimino_group = mol.HasSubstructMatch(Chem.MolFromSmarts("[O-]N=C"))

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, "Molecular weight outside typical cephalosporin range"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Insufficient rotatable bonds for cephalosporin"

    if has_amino_acid_side_chain and (has_heteroaromatic_ring or has_oxyimino_group):
        return True, "Contains beta-lactam, dihydrothiazine ring, amino acid side chain, and common substituents"

    return True, "Contains beta-lactam and dihydrothiazine ring with appropriate molecular properties"