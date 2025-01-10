"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: lipopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of a lipid A moiety (a disaccharide of glucosamine with attached fatty acids),
    core oligosaccharide containing heptose and Kdo sugars, and O-antigen polysaccharide chains.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize molecule
    try:
        Chem.SanitizeMol(mol)
    except:
        return False, "Molecule could not be sanitized"

    # Define patterns

    # Glucosamine unit (aminosugar)
    glucosamine_pattern = Chem.MolFromSmarts("NC[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O")
    # Disaccharide of glucosamine (Lipid A backbone)
    lipid_a_pattern = Chem.MolFromSmarts("NC[C@H]1O[C@H](CO)[C@@H](O)[C@H](O[C@H]2[C@@H](O)[C@@H](O)[C@H](O)[C@H](CO)O2)[C@@H]1O")
    # Fatty acid chain attached via amide linkage to glucosamine
    fatty_acid_amide_pattern = Chem.MolFromSmarts("NC(=O)CCCCC")
    # Kdo (3-deoxy-D-manno-oct-2-ulosonic acid)
    kdo_pattern = Chem.MolFromSmarts("O=C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H](COC=O)O1")
    # Heptose unit (7-membered sugar ring)
    heptose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O1")

    # Check for Lipid A backbone (disaccharide of glucosamine)
    if not mol.HasSubstructMatch(lipid_a_pattern):
        return False, "Lipid A backbone not found"

    # Check for at least four fatty acid chains attached via amide or ester linkages
    fatty_acid_ester_pattern = Chem.MolFromSmarts("C(=O)O[C;H2]")
    fatty_acid_amide_pattern = Chem.MolFromSmarts("C(=O)N[C;H2]")
    fatty_acid_chains = mol.GetSubstructMatches(fatty_acid_ester_pattern) + mol.GetSubstructMatches(fatty_acid_amide_pattern)
    if len(fatty_acid_chains) < 4:
        return False, f"Found {len(fatty_acid_chains)} fatty acid chains, need at least 4"

    # Check for Kdo units
    kdo_matches = mol.GetSubstructMatches(kdo_pattern)
    if len(kdo_matches) == 0:
        return False, "No Kdo units found"

    # Check for at least two heptose units
    heptose_matches = mol.GetSubstructMatches(heptose_pattern)
    if len(heptose_matches) < 2:
        return False, f"Found {len(heptose_matches)} heptose units, need at least 2"

    # Optional: Check for high molecular weight
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 1000:
        return False, f"Molecular weight ({mol_weight:.2f} Da) too low for lipopolysaccharide"

    return True, "Contains Lipid A backbone with attached fatty acids, Kdo, and heptose units"

__metadata__ = {
    'chemical_class': {
        'name': 'lipopolysaccharide',
        'definition': 'Liposaccharide natural compounds consisting of a trisaccharide repeating unit (two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic acid units (they are a major constituent of the cell walls of Gram-negative bacteria).'
    }
}