"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine contains an alkyl group attached to an aromatic ring and an amine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for aromatic rings
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if not aromatic_atoms:
        return None, "No aromatic ring found"

    # Check for amine groups (search for nitrogen atoms)
    nitrogen_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return None, "No nitrogen atom found, aminic group missing"

    # Check for alkyl chains
    # We look for aliphatic carbon chains
    alkyl_chain_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic()]
    if not alkyl_chain_atoms:
        return None, "No alkyl chain found"

    # Verify connectivity: check if there is a path (bond) from aromatic atoms to alkyl to nitrogen
    for ar_atom in aromatic_atoms:
        for al_atom in alkyl_chain_atoms:
            if mol.HasPath(Chem.MolFragment(mol, [ar_atom, al_atom])):
                for n_atom in nitrogen_atoms:
                    if mol.HasPath(Chem.MolFragment(mol, [al_atom, n_atom])):
                        return True, "Contains aromatic ring, alkyl chain, and amine group attached appropriately"

    return None, "Cannot verify structural features of aralkylamine"