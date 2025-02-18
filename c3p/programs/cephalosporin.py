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
    six-membered dihydrothiazine (or dihydrooxazine) ring.

    This function looks for:
    - A beta-lactam ring (four-membered cyclic amide with one nitrogen and one carbonyl group)
    - Fused to a six-membered ring containing one nitrogen and one sulfur or oxygen atom
    - The fused ring system should have specific connectivity characteristic of cephalosporins

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

    # Define SMARTS pattern for cephalosporin core
    ceph_core_smarts = """
    [#6]-1=[#6]-[#7]-2-[#6]-[#16,#8]-[#6]-[#6]-2-[#6]-1
    """

    ceph_core = Chem.MolFromSmarts(ceph_core_smarts)
    if ceph_core is None:
        return False, "Invalid SMARTS pattern for cephalosporin core"

    # Check for cephalosporin core substructure match
    if mol.HasSubstructMatch(ceph_core):
        return True, "Contains cephalosporin core structure with fused beta-lactam and dihydrothiazine/dihydrooxazine rings"
    else:
        return False, "Cephalosporin core structure not found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:23066',
        'name': 'cephalosporin',
        'definition': 'A class of beta-lactam antibiotics differing from the penicillins in having a 6-membered, rather than a 5-membered, side ring. Although cephalosporins are among the most commonly used antibiotics in the treatment of routine infections, and their use is increasing over time, they can cause a range of hypersensitivity reactions, from mild, delayed-onset cutaneous reactions to life-threatening anaphylaxis in patients with immunoglobulin E (IgE)-mediated allergy.',
    },
}