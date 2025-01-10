"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax typically contains long alkyl chains connected via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester group presence: O=C-O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkages found"

    # Identify long alkyl chains
    # This pattern looks for a chain of at least 10 carbon atoms
    long_alkyl_chain = Chem.MolFromSmarts("C(-*){10,}")
    alkyl_chain_matches = mol.GetSubstructMatches(long_alkyl_chain)
    if len(alkyl_chain_matches) < 1:
        return False, "No long alkyl chains found, required at least 10 carbon atoms"

    # Check for molecular weight to confirm large size
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Waxes are generally larger molecules
        return False, "Molecular weight too low for a typical wax compound"

    return True, "Contains characteristic ester linkage and long alkyl chains typical of waxes"