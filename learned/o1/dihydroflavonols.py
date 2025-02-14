"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is defined as any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the dihydroflavonol core SMARTS pattern without stereochemistry
    # This pattern represents the dihydroflavonol skeleton with a hydroxy group at position 3
    dihydroflavonol_smarts = 'O=C1CC(O)Oc2ccccc12'
    dihydroflavonol_pattern = Chem.MolFromSmarts(dihydroflavonol_smarts)
    if dihydroflavonol_pattern is None:
        return False, "Invalid dihydroflavonol SMARTS pattern"

    # Check if the molecule contains the dihydroflavonol core
    matches = mol.GetSubstructMatches(dihydroflavonol_pattern)
    if not matches:
        return False, "Molecule does not contain the dihydroflavonol core structure"

    # For each match, verify the hydroxy group at position 3
    for match in matches:
        # Atom indices in the SMARTS pattern:
        # 0: O (ketone oxygen at C4)
        # 1: C1 (carbonyl carbon at C4)
        # 2: C2 (C3 carbon with OH)
        # 3: C3 (C2 carbon connected to chroman oxygen)
        # 4: O (chroman oxygen)
        # 5-10: aromatic ring fused to chroman ring

        c3_idx = match[2]  # Index of the carbon at position 3 (with OH)
        atom_c3 = mol.GetAtomWithIdx(c3_idx)

        # Check if carbon at position 3 is attached to a hydroxy group
        has_oh = False
        for neighbor in atom_c3.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                if neighbor.GetTotalDegree() == 1:  # Should be -OH group
                    has_oh = True
                    break
        if has_oh:
            return True, "Molecule is a dihydroflavonol"

    return False, "Hydroxy group at position 3 not found"