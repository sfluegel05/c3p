"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is defined as a naturally occurring, basic nitrogen compound (mostly heterocyclic).
    Amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars, and antibiotics are
    not normally regarded as alkaloids. Compounds where nitrogen is exocyclic are usually classed
    as amines rather than alkaloids.

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

    # Check if any nitrogen atoms are in a ring (heterocyclic nitrogen)
    n_ring_atoms = [atom for atom in n_atoms if atom.IsInRing()]
    if not n_ring_atoms:
        return False, "No nitrogen atoms in ring (heterocycle)"

    # Exclude compounds where all nitrogen atoms are exocyclic
    if len(n_ring_atoms) == 0:
        return False, "All nitrogen atoms are exocyclic"

    # Check for basic nitrogen atoms
    basic_nitrogens = []
    for atom in n_ring_atoms:
        # Exclude quaternary nitrogen (N+)
        if atom.GetFormalCharge() != 0:
            continue
        # Pyrrole-like nitrogens are not basic (lone pair involved in aromaticity)
        # Pyridine-like nitrogens are basic
        # We'll consider sp2 nitrogens in aromatic rings not basic
        if atom.GetIsAromatic() and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            # Count the number of pi electrons contributed by nitrogen
            num_pi_electrons = rdMolDescriptors.CalcAtomElectronContributionToAromaticRing(atom)
            if num_pi_electrons == 1:  # Pyrrole-like nitrogen contributes lone pair to aromaticity
                continue  # Not basic
            else:
                basic_nitrogens.append(atom)
        else:
            # Nitrogen is not aromatic or not sp2 hybridized, considered basic
            basic_nitrogens.append(atom)

    if not basic_nitrogens:
        return False, "No basic nitrogen atoms found"

    # Exclude amino acids and peptides by detecting peptide bonds
    peptide_bond = Chem.MolFromSmarts("N[C;D2](=O)C")  # Generic peptide bond pattern
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Contains peptide bond, may be a peptide or protein"

    # Exclude nucleotides and nucleic acids by detecting nucleobases
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
    amino_sugar = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1[NH2]")  # Simple amino sugar pattern
    if mol.HasSubstructMatch(amino_sugar):
        return False, "Contains amino sugar moiety"

    # Exclude antibiotics - too broad to pattern match, often peptides or complex structures
    # For simplicity, we will not attempt to exclude antibiotics here

    return True, "Molecule is classified as an alkaloid"