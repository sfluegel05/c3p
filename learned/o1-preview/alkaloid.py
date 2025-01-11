"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is a naturally occurring nitrogen-containing compound (mostly heterocyclic),
    occurring mostly in plants, but also found in bacteria, fungi, and animals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found"

    # Check for nitrogen atoms in rings (heterocyclic nitrogen)
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    has_heterocyclic_nitrogen = False
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_atomic_nums = set(atom.GetAtomicNum() for atom in ring_atoms)
        if 7 in ring_atomic_nums:
            has_heterocyclic_nitrogen = True
            break

    if not has_heterocyclic_nitrogen:
        return False, "No nitrogen atoms found in rings (no heterocyclic nitrogen)"

    # Exclude amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars, and antibiotics
    # Check for standard amino acid backbone
    amino_acid = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")
    if mol.HasSubstructMatch(amino_acid):
        return False, "Structure matches amino acid backbone"

    # Check for peptide bond pattern (N-C(=O)-C)
    peptide_bond = Chem.MolFromSmarts("N-C(=O)-C")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond)
    if len(peptide_bonds) > 3:  # Arbitrary threshold to exclude peptides/proteins
        return False, "Contains multiple peptide bonds (possible peptide or protein)"

    # Check for nucleic acid bases
    purine = Chem.MolFromSmarts("c1ncnc2ncnc12")
    pyrimidine = Chem.MolFromSmarts("c1ccncn1")
    if mol.HasSubstructMatch(purine) or mol.HasSubstructMatch(pyrimidine):
        return False, "Contains purine or pyrimidine base (possible nucleotide or nucleic acid)"

    # Check for sugar moieties (e.g., furanose or pyranose rings)
    sugar_pattern = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C1O")
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Contains sugar moiety (possible amino sugar or glycoside)"

    # Check for antibiotics common substructures (e.g., beta-lactam)
    beta_lactam = Chem.MolFromSmarts("C1C(=O)NC1")
    if mol.HasSubstructMatch(beta_lactam):
        return False, "Contains beta-lactam ring (possible antibiotic)"

    # If all checks pass, classify as alkaloid
    return True, "Molecule contains heterocyclic nitrogen characteristic of alkaloids"

__metadata__ = {   'chemical_class': {   'name': 'alkaloid',
                              'definition': 'Any of the naturally occurring, basic nitrogen compounds (mostly heterocyclic) occurring mostly in the plant kingdom, but also found in bacteria, fungi, and animals. By extension, certain neutral compounds biogenetically related to basic alkaloids are also classed as alkaloids. Amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars and antibiotics are not normally regarded as alkaloids. Compounds in which the nitrogen is  exocyclic (dopamine, mescaline, serotonin, etc.) are usually classed as amines rather than alkaloids.',
                              'parents': []},
    'message': None,
    'success': True}