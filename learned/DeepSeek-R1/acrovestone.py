"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone (a polyphenol from Acronychia pedunculata with antioxidant/antityrosinase activity)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is acrovestone based on its SMILES string.
    Acrovestone is characterized as an isoflavone derivative with glycosylation and phenolic groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if matches acrovestone characteristics, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved isoflavone core detection (3-phenylchromen-4-one structure)
    # Matches benzopyran-4-one with substituent at C3 position
    isoflavone_pattern = Chem.MolFromSmarts("O=C1C2=C(C=CC=C2)OC3=C1C=CC=C3")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core detected"

    # Detect glycosylation (O-linked to carbohydrate-like structure)
    # Looks for oxygen connected to a ring atom (sugar moiety)
    glycoside_pattern = Chem.MolFromSmarts("[O;!H0]-[C;r5,r6]")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosylation detected"

    # Count phenolic -OH groups (attached to aromatic carbons)
    phenolic_oh = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
            # Check if attached to aromatic carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    phenolic_oh += 1
                    break
    if phenolic_oh < 1:
        return False, f"Insufficient phenolic OH groups ({phenolic_oh})"

    # Verify molecular weight range based on examples (400-1000 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if not (400 < mol_wt < 1000):
        return False, f"Molecular weight {mol_wt:.1f} outside typical range"

    # Additional check for typical isoflavone oxidation pattern
    ketone_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=O")))
    if ketone_count < 1:
        return False, "Missing characteristic ketone group"

    return True, "Isoflavone core with glycosylation and phenolic groups"