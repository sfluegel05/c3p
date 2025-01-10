"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: Peptide Antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolOps

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic typically contains multiple amide groups, possibly cyclic structures,
    and may contain unusual amino acids or heterocyclic groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a potential peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count amide bonds (R-C(=O)-NR)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 3:
        return False, "Insufficient amide (peptide) bonds"

    # Check for cyclic structures - ensure there's at least one ring
    num_rings = rdMolOps.GetSSSR(mol)
    if num_rings < 1:
        return False, "No cyclic structures detected"

    # Optionally count additional features, such as hydroxyls, indicative of unusual residues
    hydroxyl_count = Descriptors.NumHDonors(mol)
    if hydroxyl_count < 1:
        return False, "No hydroxyl groups detected; unusual for peptide antibiotics"

    # Molecular weight check - peptide antibiotics are generally >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for complex peptide antibiotic"

    return True, "Structural features typical of a peptide antibiotic found"