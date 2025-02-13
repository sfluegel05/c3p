"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:36357 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide consisting of three amino acid residues
    connected by peptide bonds.

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

    # Define SMARTS pattern for tripeptide
    tripeptide_pattern = Chem.MolFromSmarts("[N;X3;H1,H2][C;X4;H3](=[O;X1])[C;X4;H1,H2,H3][N;X3;H1,H2][C;X4;H3](=[O;X1])[C;X4;H1,H2,H3][N;X3;H1,H2][C;X4;H3](=[O;X1])[C;X4;H1,H2,H3]")

    # Check for tripeptide pattern match
    if not mol.HasSubstructMatch(tripeptide_pattern):
        return False, "Does not match tripeptide substructure pattern"

    # Check molecular weight range (typically 300-600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, "Molecular weight outside typical range for tripeptides"

    # Check for reasonable number of H-bond donors/acceptors
    n_hbond_donors = rdMolDescriptors.CalcNumLipinskiHBDon(mol)
    n_hbond_acceptors = rdMolDescriptors.CalcNumLipinskiHBAcc(mol)
    if n_hbond_donors < 3 or n_hbond_donors > 9 or n_hbond_acceptors < 6 or n_hbond_acceptors > 12:
        return False, "Number of H-bond donors/acceptors outside typical range for tripeptides"

    return True, "Contains three amino acid residues connected by peptide bonds in a linear sequence"