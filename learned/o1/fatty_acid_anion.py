"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: fatty acid anion
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group.

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

    # Sanitize molecule
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.MolSanitizeException as e:
        return False, f"Sanitization failed: {str(e)}"

    # Calculate total formal charge
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if charge >= 0:
        return False, "Molecule does not have a net negative charge"

    # Look for carboxylate group (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 3:
        return False, f"Molecule has {num_carbons} carbon atoms, less than minimum required for a fatty acid"

    return True, "Molecule is a fatty acid anion with a carboxylate group and sufficient carbon chain"