"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    Guaiacols are phenols with an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a custom SMARTS pattern for ortho-methoxyphenol
    guaiacol_pattern = Chem.MolFromSmarts("c1cc(OC)c(O)c1")

    # Check each aromatic ring in the molecule
    aromatic_rings = mol.GetAromaticRings()
    for ring in aromatic_rings:
        ring_submol = Chem.PathToSubmol(mol, ring)
        if ring_submol.HasSubstructMatch(guaiacol_pattern):
            return True, "Contains ortho-methoxy group adjacent to phenol on an aromatic ring"

    # Additional atom-level checks to confirm ortho relationship
    for cycle in aromatic_rings:
        ortho_hydroxy = None
        ortho_methoxy = None
        for idx in cycle:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "O":
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "C" and any(n.GetSymbol() == "O" for n in neighbor.GetNeighbors() if n != atom):
                        ortho_hydroxy = True
                    if neighbor.GetSymbol() == "C" and any(n.GetSymbol() == "C" for n in neighbor.GetNeighbors()):
                        potential_methoxy = [n for n in neighbor.GetNeighbors() if n.GetSymbol() == "O" and any(c.GetSymbol() == "C" for c in n.GetNeighbors() if c != neighbor)]
                        if potential_methoxy:
                            ortho_methoxy = True
            if ortho_hydroxy and ortho_methoxy:
                return True, "Contains ortho-methoxy group adjacent to phenol on an aromatic ring"

    return False, "The structure does not match the guaiacol pattern"