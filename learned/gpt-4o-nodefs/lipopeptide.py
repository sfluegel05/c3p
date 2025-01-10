"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a lipopeptide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of long aliphatic chain as part of lipid moiety
    num_chains = len(rdmolops.FindAllPathsOfLengthN(mol, 16, useBonds=True))  # Example path length for fatty acid chain
    if num_chains < 1:
        return False, "No long hydrocarbon chain found (indicating lipid motif)"

    # Check for presence of amide bonds typically found in peptide chains
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:
        return False, f"Found {len(amide_matches)} amide bonds, need at least 2 for peptide linkage"

    # Count chiral centers as a proxy for peptide sequence complexity
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 3:
        return False, f"Found {len(chiral_centers)} chiral centers, indicating insufficient peptide component"

    return True, "Contains long hydrocarbon chain and multiple amide bonds and chiral centers indicative of lipopeptide structure"