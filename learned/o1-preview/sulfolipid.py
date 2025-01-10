"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonic acid group connected via sulfur-carbon bond
    # SMARTS pattern for C-S(=O)(=O)[O-] or C-S(=O)(=O)O
    # Matches both protonated and deprotonated sulfonic acid groups
    sulfonic_acid_pattern = Chem.MolFromSmarts("[#6]-[S](=O)(=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group connected to carbon found"

    # Count the number of carbon atoms to assess lipid characteristic
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms for a lipid"

    # Count rotatable bonds to assess chain length and flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too few rotatable bonds for a lipid"

    # Check molecular weight - sulfolipids typically have high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a sulfolipid"

    return True, "Contains sulfonic acid group connected via carbon-sulfur bond and lipid characteristics"