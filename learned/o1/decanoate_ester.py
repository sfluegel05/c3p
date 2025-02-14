"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation
    of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define decanoate ester SMARTS pattern
    # This pattern matches an ester group with a linear 10-carbon acyl chain (decanoyl group)
    decanoate_smarts = '[CH3][CH2]{8}C(=O)O[!$(*C(=O)[O,N])]'

    decanoate_mol = Chem.MolFromSmarts(decanoate_smarts)
    if decanoate_mol is None:
        return False, "Invalid decanoate SMARTS pattern"

    # Search for the decanoate ester pattern in the molecule
    if not mol.HasSubstructMatch(decanoate_mol):
        return False, "No decanoate ester groups found"

    # Find all matches
    matches = mol.GetSubstructMatches(decanoate_mol)
    for match in matches:
        # Extract atoms in the acyl chain
        acyl_chain_atoms = match[:-2]  # Exclude the ester oxygen and attached atom
        acyl_chain = [mol.GetAtomWithIdx(idx) for idx in acyl_chain_atoms]

        # Check that the acyl chain is unbranched and saturated
        for atom in acyl_chain:
            if atom.GetDegree() != 2:
                # Terminal carbons can have degree 1 (methyl group)
                if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1:
                    continue
                else:
                    break
            # Check for unsaturation
            bonds = [bond for bond in atom.GetBonds() if bond.GetBondType() != Chem.rdchem.BondType.SINGLE]
            if bonds:
                break
        else:
            return True, "Contains decanoate ester group"

    return False, "No decanoate ester groups found"