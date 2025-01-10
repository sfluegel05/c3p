"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group (-C([O-])=O)
    and a hydrocarbon chain, which may include double bonds or functional groups,
    and should have sufficient chain length.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Require a carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group (-C([O-])=O) found"

    # Detecting hydrocarbons: exclude aromatic/cyclic
    non_aromat_hydrocarbon_pattern = Chem.MolFromSmarts("[#6;!R]") # carbon atoms, non-aromatic, avoids cycles
    non_aromat_hydrocarbons = mol.GetSubstructMatches(non_aromat_hydrocarbon_pattern)

    # Check if there is a continuous path of at least 6 carbon atoms excluding rings
    if len(non_aromat_hydrocarbons) < 6:
        return False, "Insufficient linear carbon chain length for fatty acid anion"

    # Check for the absence of overly complex functionalization
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings > 0:
        return False, "Presence of rings, not typical for classic fatty acid anions"

    return True, "Contains a deprotonated carboxylate group with a sufficiently long carbon chain and no ring structures"