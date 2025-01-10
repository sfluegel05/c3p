"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide comprises a peptide moiety and an attached lipid (flexible long hydrocarbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify peptide bonds (amide linkages)
    peptide_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) < 2:  # Require at least 2 peptide bonds for a protein-like fragment
        return False, "Insufficient peptide bonds found"
    
    # Identify lipid-like long hydrocarbon chains
    # Broad pattern for long aliphatic chains, considering branching and unsaturation
    lipid_pattern = Chem.MolFromSmarts("CCCCCCCCC")  # Linear chain with minimum length
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No sufficient long hydrocarbon chains found"
    
    # Check molecular weight - often indicative of complex lipopeptides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight may be too low for a lipopeptide"
    
    return True, "Contains both peptide bonds and sufficient long hydrocarbon chains, characteristic of lipopeptides"