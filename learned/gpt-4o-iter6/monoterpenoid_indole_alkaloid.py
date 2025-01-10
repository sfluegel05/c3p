"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Monoterpenoid indole alkaloids are complex structures derived from tryptophan and specialized 
    diisoprenoid precursors, with characteristic indole and expanded heterocyclic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded indole-like structures covering diversity found in the class
    indole_patterns = [
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3ccccc23"),  # Standard indole
        Chem.MolFromSmarts("C1=CC2=C3C=CC=CC3=NC2=C1"),  # Fused indoles
        Chem.MolFromSmarts("C1=NC2=CC=CC=C2C1C3=CC=CC=C3"),  # Indole and analogs
        Chem.MolFromSmarts("c1cc2cnc(c2cc1)N"),  # Indolinamine
        Chem.MolFromSmarts("c1cc2[nH]cnc2c1"),  # Pyrroloindoles
        # More sophisticated and specific indole-like patterns can be included here.
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in indole_patterns):
        return False, "No recognizable indole-like structure"

    # Check for the presence of multiple nitrogen atoms, since alkaloids often feature such.
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Insufficient nitrogen atoms for typical alkaloid complexity"

    # Total ring count and aromatic complexity - expecting rich chemistry in these compounds
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_rings < 3 or num_aromatic_rings < 2:
        return False, f"Not enough total or aromatic rings: {num_rings} total, {num_aromatic_rings} aromatic"

    # Recognize isoprene unit derivatives typical in monoterpenoids
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "Missing diisoprenoid units typical of monoterpenoids"

    # Expanded recognition of functional groups relevant to the class, e.g., methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("C-OC")  # General methoxy group
    if not mol.HasSubstructMatch(methoxy_pattern):
        return False, "Missing methoxy groups characteristic to these structures"

    return True, "Contains indole core, sufficient nitrogen, terpenoid features and methoxy groups"