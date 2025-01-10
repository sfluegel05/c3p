"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: A compound in which a sulfanyl group, -SH, is attached to an alkyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Look for -SH (sulfanyl) group pattern
    sulfanyl_pattern = Chem.MolFromSmarts("[SH1]")
    if not mol.HasSubstructMatch(sulfanyl_pattern):
        return False, "No sulfanyl (-SH) group found"

    # Check if sulfur is connected to carbon
    for s_atom in sulfur_atoms:
        neighbors = s_atom.GetNeighbors()
        if len(neighbors) == 1:  # Terminal sulfur (as in -SH)
            neighbor = neighbors[0]
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if the carbon is part of an alkyl group
                # (not aromatic, not part of C=O, etc)
                if (not neighbor.GetIsAromatic() and 
                    neighbor.GetHybridization() == Chem.HybridizationType.SP3):
                    return True, "Contains sulfanyl group (-SH) attached to an alkyl group"

    return False, "Sulfanyl group not attached to an alkyl carbon"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        "SCCCC",  # butanethiol
        "CCS",    # ethanethiol
        "SC(CCC)C",  # 2-pentanethiol
        "SC1CCCC1",  # cyclopentanethiol
        "SCCCCCCC",  # heptane-1-thiol
        "CC(C)S",    # propane-2-thiol
        "C=CCS",     # 2-propene-1-thiol
        "c1ccccc1S",  # thiophenol (not an alkanethiol)
        "CCSSC",      # diethyl disulfide (not an alkanethiol)
        "CC(=O)S"     # thioacetic acid (not an alkanethiol)
    ]
    
    for smiles in examples:
        result, reason = is_alkanethiol(smiles)
        print(f"SMILES: {smiles}")
        print(f"Is alkanethiol: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()