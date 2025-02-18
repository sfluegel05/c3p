"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: CHEBI:23066 cephalosporin
"""

from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are beta-lactam antibiotics characterized by a beta-lactam ring fused to a 
    dihydrothiazine ring (a six-membered ring containing one sulfur and one nitrogen atom).
    
    This function looks for:
    - A beta-lactam ring: a four-membered cyclic amide (lactam) ring with one nitrogen atom.
    - A dihydrothiazine ring: a six-membered ring with one nitrogen and one sulfur atom.
    - The beta-lactam ring and dihydrothiazine ring are fused (share two atoms and a bond).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-lactam ring (4-membered cyclic amide with one nitrogen atom)
    beta_lactam_smarts = '[R2,N;R][C;R]=[O;R,R1][C;R]'

    beta_lactam = Chem.MolFromSmarts(beta_lactam_smarts)
    if beta_lactam is None:
        return False, "Invalid SMARTS pattern for beta-lactam ring"

    # Define SMARTS pattern for dihydrothiazine ring (6-membered ring with one nitrogen and one sulfur atom)
    dihydrothiazine_smarts = 'C1C(=C)SC(N)C=C1'

    dihydrothiazine = Chem.MolFromSmarts(dihydrothiazine_smarts)
    if dihydrothiazine is None:
        return False, "Invalid SMARTS pattern for dihydrothiazine ring"

    # Find beta-lactam matches
    beta_lactam_matches = mol.GetSubstructMatches(beta_lactam)
    if not beta_lactam_matches:
        return False, "Beta-lactam ring not found"

    # Find dihydrothiazine matches
    dihydrothiazine_matches = mol.GetSubstructMatches(dihydrothiazine)
    if not dihydrothiazine_matches:
        return False, "Dihydrothiazine ring not found"

    # Check if beta-lactam and dihydrothiazine rings are fused
    fused = False
    for bl_ring in beta_lactam_matches:
        bl_atoms = set(bl_ring)
        for dt_ring in dihydrothiazine_matches:
            dt_atoms = set(dt_ring)
            shared_atoms = bl_atoms & dt_atoms
            if len(shared_atoms) >= 2:
                # Rings are fused if they share at least two atoms
                fused = True
                break
        if fused:
            break

    if not fused:
        return False, "Beta-lactam ring and dihydrothiazine ring are not fused"

    return True, "Contains cephalosporin core structure with fused beta-lactam and dihydrothiazine rings"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:23066',
        'name': 'cephalosporin',
        'definition': 'A class of beta-lactam antibiotics differing from the penicillins in having a 6-membered, rather than a 5-membered, side ring. Although cephalosporins are among the most commonly used antibiotics in the treatment of routine infections, and their use is increasing over time, they can cause a range of hypersensitivity reactions, from mild, delayed-onset cutaneous reactions to life-threatening anaphylaxis in patients with immunoglobulin E (IgE)-mediated allergy.',
    },
}