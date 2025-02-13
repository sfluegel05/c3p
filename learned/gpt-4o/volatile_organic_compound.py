"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is assumed to have an initial boiling point <= 250°C under standard atmospheric pressure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a VOC, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)

    # Count number of carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for specific functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    alkene_pattern = Chem.MolFromSmarts("C=C")
    alkyne_pattern = Chem.MolFromSmarts("C#C")
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I,F]")

    # Evaluate presence and potential impact of functional groups
    volatility_based_on_groups = (
        mol.HasSubstructMatch(alcohol_pattern) or
        mol.HasSubstructMatch(ether_pattern) or
        mol.HasSubstructMatch(alkene_pattern) or
        mol.HasSubstructMatch(alkyne_pattern) or
        mol.HasSubstructMatch(halogen_pattern)
    )

    # Use adjusted criteria for a more accurate classification
    # Focus on functional groups and not exceeding a high molecular weight threshold
    if volatility_based_on_groups and mol_wt <= 350:  # Adjusted MW threshold to align closer with expected VOC volatility
        return True, "Contains VOC functional group and structure suggests potential volatility"

    # Consider molecules with up to 20 carbons as possibly VOCs, barring complexity from heavy atoms or rings
    if carbon_count <= 20:
        return True, "Carbon chain length and structure suggest potential volatility"

    return False, "Structure does not indicate VOC features"