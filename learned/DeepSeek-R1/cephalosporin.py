"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: Cephalosporin antibiotics (CHEBI:23066)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins contain a beta-lactam ring fused to a dihydrothiazine ring (6-membered sulfur-containing ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Core bicyclo[4.2.0]oct-2-ene system with sulfur and beta-lactam
    core_pattern = Chem.MolFromSmarts("[C]1([SX2][CX4][CX4][CX4][CX4]1)[NX3][C](=O)[C]2=C([CX4])[C](=O)N12")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing bicyclic beta-lactam core"
    
    # Verify beta-lactam (4-membered ring with amide)
    beta_lactam = Chem.MolFromSmarts("[NX3][C](=O)[CX4][CX4]1")
    beta_lactam_matches = mol.GetSubstructMatches(beta_lactam)
    if not beta_lactam_matches:
        return False, "No beta-lactam ring found"
    
    # Check for sulfur in the six-membered ring
    sulfur_in_ring = False
    s_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[S]"))
    for s_idx in [a[0] for a in s_matches]:
        atom = mol.GetAtomWithIdx(s_idx)
        if atom.IsInRingSize(6):
            sulfur_in_ring = True
            break
    if not sulfur_in_ring:
        return False, "Sulfur not in six-membered ring"
    
    return True, "Contains beta-lactam fused to dihydrothiazine ring"