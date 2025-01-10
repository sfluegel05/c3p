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
    # Let's check for paths of length 14 to 28, representing longer fatty acids
    long_chain_found = False
    for length in range(14, 29):
        num_chains = len(rdmolops.FindAllPathsOfLengthN(mol, length, useBonds=True))
        if num_chains > 0:
            long_chain_found = True
            break

    if not long_chain_found:
        return False, "No long hydrocarbon chain found (indicating lipid motif)"

    # Check for presence of amide bonds typically found in peptide chains
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 4:  # Increased to a minimum of 4
        return False, f"Found {len(amide_matches)} amide bonds, need at least 4 for peptide linkage"
    
    # Count chiral centers as a proxy for peptide sequence complexity
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 3:
        return False, f"Found {len(chiral_centers)} chiral centers, indicating insufficient peptide component"
    
    # Attempt to detect an amphiphilic pattern:
    # Check if there are at least 2 distinct hydrophobic and hydrophilic zones
    num_hydrophobic = len(list(filter(lambda x: x.GetAtomicNum() == 6, mol.GetAtoms())))
    num_hydrophilic = len(list(filter(lambda x: x.GetAtomicNum() in [7, 8, 16], mol.GetAtoms())))  # N and O primarily
    
    if num_hydrophobic < 12 or num_hydrophilic < 3:
        return False, "Does not have enough hydrophobic and hydrophilic regions for a typical amphiphilic lipopeptide"

    return True, "Contains amphiphilic character with long hydrocarbon chain, multiple amide bonds, and several chiral centers indicative of lipopeptide structure"