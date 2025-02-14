"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:35346 sulfolipid
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

    # Define a SMARTS pattern for sulfonic acid residue connected via C-S bond
    # Sulfonic acid group: S(=O)(=O)-O
    # Connected via sulfur to a carbon atom
    sulfonic_acid_pattern = Chem.MolFromSmarts("[#6]-S(=O)(=O)-[O;H1]")
    if sulfonic_acid_pattern is None:
        return False, "Invalid SMARTS pattern for sulfonic acid residue"

    # Check for the sulfonic acid group connected via C-S bond
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid residue connected via carbon-sulfur bond found"

    # Optionally, check for lipid-like characteristics (long hydrocarbon chains)
    # Count the number of carbons in aliphatic chains (excluding rings and aromatic carbons)
    aliphatic_carbons = [atom for atom in mol.GetAtoms() 
                         if atom.GetAtomicNum() == 6 and not atom.IsInRing() and not atom.GetIsAromatic()]
    if len(aliphatic_carbons) < 10:
        return False, "Too few aliphatic carbons for a lipid"

    # Check for the presence of long carbon chains (length >= 10)
    chains = Chem.rdmolops.GetLongestAliphaticChain(mol)
    if chains < 10:
        return False, "No long hydrocarbon chains found"

    return True, "Contains sulfonic acid residue connected via carbon-sulfur bond to a lipid"