"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester

Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is an ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester pattern with labeled atoms for acyl carbon (atom 1) and ester oxygen (atom 2)
    ester_pattern = Chem.MolFromSmarts('[C:1](=O)[O:2][C]')

    # Find ester matches
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Define butyryl group (butanoic acid minus the OH)
    butyryl_smiles = 'CCCC(=O)'
    butyryl_mol = Chem.MolFromSmiles(butyryl_smiles)

    # For each ester, extract the acyl fragment and check if it matches butyryl
    for match in ester_matches:
        acyl_carbon_idx = match[0]  # Atom index of the acyl carbon (C=O)
        ester_oxygen_idx = match[1]  # Atom index of the ester oxygen

        # Create a copy of the molecule to modify
        mol_copy = Chem.RWMol(mol)

        # Remove the bond between acyl carbon and ester oxygen to isolate acyl part
        mol_copy.RemoveBond(acyl_carbon_idx, ester_oxygen_idx)

        # Get the fragments resulting from bond cleavage
        frags = Chem.GetMolFrags(mol_copy.GetMol(), asMols=True, sanitizeFrags=True)

        # Find the fragment containing the acyl carbon
        for frag in frags:
            atom_indices = [atom.GetAtomMapNum() for atom in frag.GetAtoms()]
            # Set atom map numbers to match indices
            for atom in frag.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx())
            if any(atom.GetAtomicNum() == 6 and atom.GetIdx() == acyl_carbon_idx for atom in frag.GetAtoms()):
                # Found the acyl fragment
                # Compare to butyryl group using substructure match
                if frag.HasSubstructMatch(butyryl_mol):
                    return True, "Contains butyrate ester group"

    return False, "Does not contain butyrate ester group"