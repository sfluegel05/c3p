"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35627 beta-lactam antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

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

    # Aromatize molecule to correctly perceive rings
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # Define SMARTS patterns for different beta-lactam cores
    beta_lactam_smarts = '[N;R][C;R](=O)[C;R][C;R]'  # Beta-lactam ring
    penam_smarts = '[N;R1][C;R1](=O)[C;R1]1[C;R1][C;R1][S;R1]1'  # Penam core
    cephem_smarts = '[N;R1][C;R1](=O)[C;R1]1[C;R1]=[C;R1][C;R1][S;R1]1'  # Cephem core
    carbapenem_smarts = '[N;R1][C;R1](=O)[C;R1]1=[C;R1][C;R1][C;R1]1'  # Carbapenem core
    clavam_smarts = '[N;R1][C;R1](=O)[C;R1]1[C;R1][O;R1][C;R1]1'  # Clavam core
    monobactam_smarts = '[N;R1][C;R1](=O)[C;R1][C;R1]'  # Monobactam core (non-fused beta-lactam)

    # Compile patterns
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    penam_pattern = Chem.MolFromSmarts(penam_smarts)
    cephem_pattern = Chem.MolFromSmarts(cephem_smarts)
    carbapenem_pattern = Chem.MolFromSmarts(carbapenem_smarts)
    clavam_pattern = Chem.MolFromSmarts(clavam_smarts)
    monobactam_pattern = Chem.MolFromSmarts(monobactam_smarts)

    # Check for beta-lactam ring
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Check for specific beta-lactam antibiotic cores
    is_penam = mol.HasSubstructMatch(penam_pattern)
    is_cephem = mol.HasSubstructMatch(cephem_pattern)
    is_carbapenem = mol.HasSubstructMatch(carbapenem_pattern)
    is_clavam = mol.HasSubstructMatch(clavam_pattern)
    is_monobactam = mol.HasSubstructMatch(monobactam_pattern)
    
    if not any([is_penam, is_cephem, is_carbapenem, is_clavam, is_monobactam]):
        return False, "Beta-lactam ring not part of known antibiotic cores"

    # Check for organonitrogen heterocyclic nature
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    is_organonitrogen_heterocycle = False

    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check if ring contains only non-metal elements and at least one nitrogen
        if any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            if all(atom.GetAtomicNum() > 1 and atom.GetAtomicNum() <= 16 for atom in ring_atoms):
                is_organonitrogen_heterocycle = True
                break

    if not is_organonitrogen_heterocycle:
        return False, "Molecule is not an organonitrogen heterocycle"

    # Check for carboxylic acid group (common in beta-lactam antibiotics)
    carboxylic_acid_smarts = '[CX3](=O)[OX1H0-]'
    carboxylic_acid_pattern = Chem.MolFromSmarts(carboxylic_acid_smarts)
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains beta-lactam antibiotic core and is an organonitrogen heterocycle"