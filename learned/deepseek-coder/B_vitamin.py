"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: CHEBI:33284 B vitamin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    B vitamins are a group of water-soluble vitamins with specific structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitrogen-containing heterocycles (common in B vitamins)
    heterocycle_patterns = [
        "[n]1ccccc1",  # Pyridine
        "[n]1ccncc1",  # Pyrimidine
        "[n]1cncc1",   # Imidazole
        "[n]1ccnc1",   # Pyrazole
        "[n]1cc[nH]c1", # Pyrrole
        "[n]1c[nH]cc1", # Pyrimidine (alternative)
        "[n]1c[nH]cn1", # Imidazole (alternative)
        "[n]1cc[nH]c1", # Pyridine (alternative)
        "[n]1cc[nH]cc1", # Pyrimidine (alternative)
        "[n]1cc[nH]cn1", # Imidazole (alternative)
        "[n]1c2ccccc2cc1", # Quinoxaline
        "[n]1c2ccccc2nc1", # Quinoline
        "[n]1c2ccccc2cn1", # Isoquinoline
        "[n]1c2ccccc2n1", # Indole
        "[n]1c2ccccc2c1", # Benzimidazole
        "[n]1c2ccccc2c1", # Purine
        "[n]1c2ccccc2c1", # Pteridine
        "[n]1c2ccccc2c1", # Flavin
        "[n]1c2ccccc2c1", # Thiazole
        "[n]1c2ccccc2c1"  # Isoalloxazine
    ]
    has_heterocycle = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in heterocycle_patterns)
    if not has_heterocycle:
        return False, "No nitrogen-containing heterocycle found"

    # Check for functional groups common in B vitamins
    functional_group_patterns = [
        "[CX3](=O)[OX2H1]",  # Carboxylic acid
        "[NX3][CX4]",        # Amine
        "[PX4](=O)([OX2H1])",# Phosphate
        "[CX3](=O)[OX1H0-]", # Carboxylate
        "[CX3](=O)[NX3]",    # Amide
        "[SX2](=O)(=O)[OX2H1]", # Sulfonate
        "[OX2H1]",           # Hydroxyl
        "[NX3H2]",           # Primary amine
        "[NX3H1]",           # Secondary amine
        "[NX3H0]",           # Tertiary amine
        "[CX3](=O)[OX2H0-]", # Ester
        "[CX3](=O)[NX3H1]",  # Amide
        "[CX3](=O)[NX3H0]",  # Amide
        "[CX3](=O)[NX3H2]",  # Amide
        "[CX3](=O)[NX3H0]",  # Amide
        "[CX3](=O)[NX3H1]",  # Amide
        "[CX3](=O)[NX3H2]",  # Amide
        "[CX3](=O)[NX3H0]",  # Amide
        "[CX3](=O)[NX3H1]",  # Amide
        "[CX3](=O)[NX3H2]"   # Amide
    ]
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in functional_group_patterns)
    if not has_functional_group:
        return False, "No functional group common in B vitamins found"

    # Check for water solubility (B vitamins are water-soluble)
    logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if logP > 1.0:  # More strict threshold for water solubility
        return False, "Molecule is likely not water-soluble"

    # Check molecular weight (B vitamins typically have MW < 1500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1500:  # Relaxed threshold for molecular weight
        return False, "Molecular weight too high for a B vitamin"

    # Count polar atoms (B vitamins typically have many polar atoms)
    polar_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8])
    if polar_atom_count < 4:
        return False, "Too few polar atoms for a B vitamin"

    return True, "Contains nitrogen-containing heterocycle, functional groups common in B vitamins, and is water-soluble"