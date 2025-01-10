"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is a naturally occurring, basic nitrogen compound (mostly heterocyclic).

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
    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_nitrogen:
        return False, "No nitrogen atoms found"

    # Check for nitrogen in heterocyclic rings
    rs = mol.GetRingInfo()
    atom_rings = rs.AtomRings()
    heterocyclic_nitrogen = False
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            # Check if ring contains at least one non-nitrogen atom (heterocycle)
            if any(atom.GetAtomicNum() != 7 for atom in ring_atoms):
                heterocyclic_nitrogen = True
                break

    if not heterocyclic_nitrogen:
        return False, "No nitrogen found in heterocyclic rings"

    # Check for basic nitrogen (sp3-hybridized nitrogen capable of protonation)
    basic_nitrogen = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Exclude quaternary ammonium (already positively charged)
            if atom.GetFormalCharge() == 1:
                continue
            # Exclude amide nitrogens (connected to a carbonyl group)
            is_amide = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    for sub_neighbor in neighbor.GetNeighbors():
                        if sub_neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(neighbor.GetIdx(), sub_neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            is_amide = True
                            break
                if is_amide:
                    break
            if is_amide:
                continue
            # Check if nitrogen is sp3-hybridized
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                basic_nitrogen = True
                break

    if not basic_nitrogen:
        return False, "No basic (sp3-hybridized) nitrogen atoms found"

    # Exclude amino acids, peptides, nucleotides, amino sugars, and antibiotics
    # Check for peptide bond pattern (N-C(=O)-C)
    peptide_bond = Chem.MolFromSmarts("N-C(=O)-C")
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Contains peptide bonds (possible peptide or protein)"

    # Check for sugar moieties (e.g., furanose or pyranose rings)
    sugar_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O")
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Contains sugar moiety (possible glycoside or amino sugar)"

    # Check for nucleotide bases (purines and pyrimidines)
    purine = Chem.MolFromSmarts("c1ncnc2ncnc12")
    pyrimidine = Chem.MolFromSmarts("c1cncnc1")
    if mol.HasSubstructMatch(purine) or mol.HasSubstructMatch(pyrimidine):
        return False, "Contains purine or pyrimidine base (possible nucleotide)"

    # Check for antibiotic-like structures (complex polyketides, non-ribosomal peptides)
    # This is complex; for simplicity, we may skip or assume acceptable unless known patterns are found

    # Check if nitrogen is exocyclic (connected outside a ring)
    exocyclic_nitrogen = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            if not atom.IsInRing():
                exocyclic_nitrogen = True
                break
    if exocyclic_nitrogen and not heterocyclic_nitrogen:
        return False, "Nitrogen is exocyclic (compound may be an amine, not an alkaloid)"

    # If all checks pass, the molecule is classified as an alkaloid
    return True, "Molecule contains basic nitrogen in a heterocyclic ring characteristic of alkaloids"

__metadata__ = {   'chemical_class': {   'name': 'alkaloid',
                              'definition': 'Any of the naturally occurring, basic nitrogen compounds (mostly heterocyclic) occurring mostly in the plant kingdom, but also found in bacteria, fungi, and animals. By extension, certain neutral compounds biogenetically related to basic alkaloids are also classed as alkaloids. Amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars and antibiotics are not normally regarded as alkaloids. Compounds in which the nitrogen is exocyclic (dopamine, mescaline, serotonin, etc.) are usually classed as amines rather than alkaloids.',
                              'parents': []},
        'message': None,
        'success': True}