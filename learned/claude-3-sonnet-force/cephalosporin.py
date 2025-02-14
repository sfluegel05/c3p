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
    Cephalosporins are a class of beta-lactam antibiotics with a 6-membered side ring.

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

    # Look for fused 6-membered ring (cephalosporin pattern)
    ceph_pattern = Chem.MolFromSmarts("[C&!R]1[N,O]C2C(=O)N2C(=O)C1")
    if not mol.HasSubstructMatch(ceph_pattern):
        return False, "No 6-membered cephalosporin ring found"

    # Check for common cephalosporin substituents
    diazine_pattern = Chem.MolFromSmarts("c1nncn1")
    pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
    methoxyimino_pattern = Chem.MolFromSmarts("CON=C")
    vinyl_pattern = Chem.MolFromSmarts("C=C")

    has_diazine = mol.HasSubstructMatch(diazine_pattern)
    has_pyridine = mol.HasSubstructMatch(pyridine_pattern)
    has_methoxyimino = mol.HasSubstructMatch(methoxyimino_pattern)
    has_vinyl = mol.HasSubstructMatch(vinyl_pattern)

    if has_diazine or has_pyridine or has_methoxyimino or has_vinyl:
        return True, "Contains beta-lactam and 6-membered cephalosporin ring, with common substituents"

    # Additional checks based on molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1200:
        return False, "Molecular weight outside typical cephalosporin range"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Insufficient rotatable bonds for cephalosporin"

    return True, "Contains beta-lactam and 6-membered cephalosporin ring with appropriate molecular properties"