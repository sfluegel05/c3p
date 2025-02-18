"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:48030 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of three amino acid residues connected by two peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all amide bonds (peptide bonds)
    # SMARTS for amide bond: [CX3](=O)[NX3H0]
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H0]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 2:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 2"
    
    # Check for three alpha carbons (each connected to an amide N and a carbonyl)
    # Alpha carbon is connected to an amide nitrogen and a carbonyl group
    # SMARTS: [C]([NX3])[CX3](=O)
    alpha_carbon_pattern = Chem.MolFromSmarts("[C]([NX3])[CX3](=O)")
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_carbons) < 3:
        return False, f"Found {len(alpha_carbons)} alpha carbons, need at least 3"
    
    # Check for three amino acid residues (each with a side chain)
    # This is approximated by checking for three distinct R groups attached to alpha carbons
    # SMARTS for amino acid residue: [C]([NX3])([CX3](=O))[*]
    residue_pattern = Chem.MolFromSmarts("[C]([NX3])([CX3](=O))[!H]")
    residue_matches = mol.GetSubstructMatches(residue_pattern)
    if len(residue_matches) < 3:
        return False, f"Found {len(residue_matches)} amino acid residues, need at least 3"
    
    # Check molecular weight (tripeptides typically > 200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for tripeptide"
    
    return True, "Contains three amino acid residues connected by two amide bonds"