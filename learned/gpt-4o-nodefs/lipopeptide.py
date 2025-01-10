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
        return (False, "Invalid SMILES string")

    # Detecting lipid moiety: Long aliphatic chain and a polar group for head
    lipid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if len(lipid_matches) == 0:
        return (False, "No long hydrocarbon chain found (indicating lipid motif)")

    # Detecting multiple amide bonds, revisiting the count to be more inclusive
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 3:  # Changed to minimum 3
        return (False, f"Found {len(amide_matches)} amide bonds, need at least 3 for peptide linkage")

    # Examine chiral centers as an indication of complexity in structure
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 3:
        return (False, f"Found {len(chiral_centers)} chiral centers, indicating insufficient peptide component")
    
    # Ensure molecule has both hydrophobic segments and polar (hydrophilic) features:
    num_hydrophobic = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    num_hydrophilic = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8, 16]])  # N and O primarily
    
    if num_hydrophobic < 10 or num_hydrophilic < 3:
        return (False, "Does not have enough hydrophobic and hydrophilic regions for a typical amphiphilic lipopeptide")

    return (True, "Contains amphiphilic character with long hydrocarbon chain, multiple amide bonds, and several chiral centers indicative of lipopeptide structure")