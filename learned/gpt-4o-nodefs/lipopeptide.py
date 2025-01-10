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
    # Let's check for paths of length 12 to 24, which are representative of typical fatty acids
    long_chain_found = False
    for length in range(12, 25):
        num_chains = len(rdmolops.FindAllPathsOfLengthN(mol, length, useBonds=True))
        if num_chains > 0:
            long_chain_found = True
            break

    if not long_chain_found:
        return False, "No long hydrocarbon chain found (indicating lipid motif)"

    # Check for presence of amide bonds typically found in peptide chains
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 3:  # Increased the requirement to a minimum of 3
        return False, f"Found {len(amide_matches)} amide bonds, need at least 3 for peptide linkage"

    # Count chiral centers as a proxy for peptide sequence complexity
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 2:  # Reduced to 2, as lipopeptides may be simple in structure
        return False, f"Found {len(chiral_centers)} chiral centers, indicating insufficient peptide component"

    # This additional pattern can help identify more peptide-like properties
    # e.g., looking for typical N-terminus
    # nitrogen_pattern = Chem.MolFromSmarts("[NH2,NH3+,NH+]")

    return True, "Contains long hydrocarbon chain and multiple amide bonds and chiral centers indicative of lipopeptide structure"