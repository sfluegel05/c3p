"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is defined as a naturally occurring, basic nitrogen compound (mostly heterocyclic),
    occurring mostly in the plant kingdom, but also found in bacteria, fungi, and animals.
    By extension, certain neutral compounds biogenetically related to basic alkaloids are also classed as alkaloids.
    Amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars, and antibiotics are
    not normally regarded as alkaloids. Compounds where nitrogen is exocyclic (e.g., dopamine, mescaline,
    serotonin) are usually classed as amines rather than alkaloids.

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
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not n_atoms:
        return False, "No nitrogen atoms found"

    # Check for nitrogen atoms in rings
    n_in_ring = [atom for atom in n_atoms if atom.IsInRing()]
    if not n_in_ring:
        return False, "No nitrogen atoms in rings found, nitrogen is exocyclic"

    # Check for basic nitrogen atoms in rings
    basic_nitrogens = []
    for atom in n_in_ring:
        if atom.GetFormalCharge() != 0:
            continue  # Exclude charged nitrogens (e.g., quaternary ammonium)
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            # SP3 hybridized nitrogen in ring (e.g., piperidine) is basic
            basic_nitrogens.append(atom)
        elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            if atom.GetIsAromatic():
                num_h = atom.GetTotalNumHs()
                # Aromatic nitrogen in ring with zero hydrogens (pyridine-like) is basic
                if num_h == 0:
                    basic_nitrogens.append(atom)
            else:
                # Non-aromatic SP2 nitrogen in ring may be basic (e.g., imine nitrogen)
                basic_nitrogens.append(atom)

    if not basic_nitrogens:
        return False, "No basic nitrogen atoms in rings found"

    # Exclude amino acids and peptides by detecting peptide bonds (amide linkage)
    peptide_bond = Chem.MolFromSmarts("N[C;D2](=O)C")  # N-C(=O)-C pattern
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Contains peptide bond, may be a peptide or protein"

    # Exclude nucleotides and nucleic acids by detecting nucleobases and sugar-phosphate backbone
    nucleic_acid_bases = [
        Chem.MolFromSmarts("c1nc[nH]c(=O)[nH]1"),  # Cytosine
        Chem.MolFromSmarts("c1cc(=O)[nH]c(=O)[nH]1"),  # Uracil/Thymine
        Chem.MolFromSmarts("c1ncnc2ncnn12"),  # Adenine
        Chem.MolFromSmarts("c1[nH]c2c(n1)nc(nc2)N"),  # Guanine
    ]
    for base in nucleic_acid_bases:
        if mol.HasSubstructMatch(base):
            return False, "Contains nucleic acid base"

    # Exclude amino sugars by detecting sugar rings with amino groups
    amino_sugar = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](N)[C@@H](O)[C@H](O)[C@H]1O")
    if mol.HasSubstructMatch(amino_sugar):
        return False, "Contains amino sugar moiety"

    # Exclude antibiotics - often contain characteristic lactone or glycopeptide structures
    antibiotic_patterns = [
        Chem.MolFromSmarts("C(=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # General lactone ring
        Chem.MolFromSmarts("C1=CNC(=O)N=C1"),  # Beta-lactam ring
    ]
    for pattern in antibiotic_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains antibiotic-like structure"

    # Exclude other specific functional groups associated with excluded classes
    excluded_groups = [
        Chem.MolFromSmarts("P(=O)(O)O"),    # Phosphate group in nucleotides
        Chem.MolFromSmarts("[SX4](=O)(=O)(O)O"),  # Sulfate group
        Chem.MolFromSmarts("C(=O)N[H]"),   # Amide group (primary)
    ]
    for group in excluded_groups:
        if mol.HasSubstructMatch(group):
            return False, "Contains excluded functional group, may be an excluded compound"

    # At this point, we have a molecule with at least one basic nitrogen atom in a ring
    # that is not excluded by the other criteria
    return True, "Molecule is classified as an alkaloid"