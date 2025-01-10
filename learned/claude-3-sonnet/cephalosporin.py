"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin antibiotics
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins contain a beta-lactam ring fused to a 6-membered dihydrothiazine ring.

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

    # First check for the core cephalosporin structure:
    # Beta-lactam (4-membered ring) fused to 6-membered dihydrothiazine ring
    ceph_core = Chem.MolFromSmarts('[#16]-1-[#6]-[#6]=2-[#6](=[#6]-1)-[#7]-3-[#6](=[#8])-[#6]-[#6]-2-3')
    if not mol.HasSubstructMatch(ceph_core):
        return False, "Missing cephalosporin core structure (fused beta-lactam and dihydrothiazine rings)"

    # Check for beta-lactam ring specifically
    beta_lactam = Chem.MolFromSmarts('[#7]-1-[#6](=[#8])-[#6]-[#6]-1')
    if not mol.HasSubstructMatch(beta_lactam):
        return False, "Missing beta-lactam ring"

    # Check for carboxylic acid/carboxylate group at position 2
    carboxyl_patterns = [
        Chem.MolFromSmarts('C(=O)[OH]'),  # carboxylic acid
        Chem.MolFromSmarts('C(=O)[O-]')   # carboxylate
    ]
    has_carboxyl = any(mol.HasSubstructMatch(p) for p in carboxyl_patterns if p is not None)
    if not has_carboxyl:
        return False, "Missing carboxylic acid/carboxylate group"

    # Check for amide group at position 7
    amide = Chem.MolFromSmarts('[NH]-[#6](=[#8])')
    if not mol.HasSubstructMatch(amide):
        return False, "Missing amide group at position 7"

    # Additional checks for common substituents that might indicate a cephalosporin
    substituent_patterns = {
        'aminothiazole': Chem.MolFromSmarts('c1csc(N)n1'),
        'tetrazole': Chem.MolFromSmarts('n1nnn[nH]1'),
        'methyltetrazole': Chem.MolFromSmarts('Cn1nnnc1'),
        'oxime': Chem.MolFromSmarts('C=NOC'),
        'acetoxymethyl': Chem.MolFromSmarts('CC(=O)OC'),
        'carboxamide': Chem.MolFromSmarts('NC(=O)')
    }
    
    found_substituents = []
    for name, pattern in substituent_patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_substituents.append(name)

    # Verify ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count for cephalosporin structure"

    # Check molecular complexity
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15:
        return False, "Molecule too small to be a cephalosporin"

    # Additional check for characteristic sulfur atom
    sulfur_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16])
    if sulfur_count < 1:
        return False, "Missing characteristic sulfur atom"

    reason = "Contains cephalosporin core structure with:"
    reason += "\n- Beta-lactam fused to dihydrothiazine ring"
    reason += "\n- Carboxylic acid/carboxylate group"
    reason += "\n- Amide group"
    if found_substituents:
        reason += f"\n- Found substituents: {', '.join(found_substituents)}"

    return True, reason