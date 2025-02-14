"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35627 beta-lactam antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic is an organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for beta-lactam ring (4-membered cyclic amide)
    beta_lactam_smarts = '[NX3;H1;R1][C;R1](=O)[C;R1][C;R1]'  # Specific to beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Define SMARTS patterns for specific beta-lactam antibiotic cores
    penam_smarts = '[NX3;H1;R1][C;R1](=O)[C;R1]1[C;R1][C;R1][S;R1]1'  # Penam core
    cephem_smarts = '[NX3;H1;R1][C;R1](=O)[C;R1]1=[C;R1][C;R1][S;R1]1'  # Cephem core
    carbapenem_smarts = '[NX3;H1;R1][C;R1](=O)[C;R1]1=[C;R1][C;R1][C;R1]1'  # Carbapenem core
    oxapenam_smarts = '[NX3;H1;R1][C;R1](=O)[C;R1]1[C;R1][O;R1][C;R1]1'  # Clavams and oxapenams
    monobactam_smarts = '[NX3;H1;R1][C;R1](=O)[C;R1][C;R1]'  # Monobactams

    # Compile patterns
    penam_pattern = Chem.MolFromSmarts(penam_smarts)
    cephem_pattern = Chem.MolFromSmarts(cephem_smarts)
    carbapenem_pattern = Chem.MolFromSmarts(carbapenem_smarts)
    oxapenam_pattern = Chem.MolFromSmarts(oxapenam_smarts)
    monobactam_pattern = Chem.MolFromSmarts(monobactam_smarts)

    # Check for specific beta-lactam antibiotic cores
    is_penam = mol.HasSubstructMatch(penam_pattern)
    is_cephem = mol.HasSubstructMatch(cephem_pattern)
    is_carbapenem = mol.HasSubstructMatch(carbapenem_pattern)
    is_oxapenam = mol.HasSubstructMatch(oxapenam_pattern)
    is_monobactam = mol.HasSubstructMatch(monobactam_pattern)

    if not any([is_penam, is_cephem, is_carbapenem, is_oxapenam, is_monobactam]):
        return False, "Beta-lactam ring not part of known antibiotic cores"

    # Verify that the molecule is an organonitrogen heterocycle
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    organonitrogen_heterocycle = False

    for ring in rings:
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check if ring contains nitrogen and is all non-metal elements
        if any(atom.GetAtomicNum() == 7 for atom in atoms_in_ring):
            if all(atom.GetAtomicNum() > 0 and atom.GetAtomicNum() <= 16 for atom in atoms_in_ring):
                organonitrogen_heterocycle = True
                break

    if not organonitrogen_heterocycle:
        return False, "Molecule is not an organonitrogen heterocycle"

    # Check for carboxylic acid group or its esters (common in beta-lactam antibiotics)
    carboxylic_acid_smarts = '[CX3](=O)[OX1H0-,OX2H1,OX2H0][#0,#6]'
    carboxylic_acid_pattern = Chem.MolFromSmarts(carboxylic_acid_smarts)
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group or ester found"

    return True, "Contains beta-lactam antibiotic core and is an organonitrogen heterocycle"