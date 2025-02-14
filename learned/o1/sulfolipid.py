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
    A sulfolipid is a compound containing a sulfonic acid residue joined by a 
    carbon-sulfur bond to a lipid.

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

    # Define SMARTS pattern for sulfonic acid group connected via carbon-sulfur bond
    sulfo_pattern = Chem.MolFromSmarts("C-S(=O)(=O)-O[H1]")
    if sulfo_pattern is None:
        return False, "Invalid SMARTS pattern for sulfonic acid group"

    # Check for the sulfonic acid group connected via carbon-sulfur bond
    sulfo_matches = mol.GetSubstructMatches(sulfo_pattern)
    if not sulfo_matches:
        return False, "No sulfonic acid residue connected via carbon-sulfur bond found"

    # Check for long aliphatic carbon chains (lipid chains)
    # Define SMARTS pattern for long aliphatic chain (at least 10 carbons)
    lipid_pattern = Chem.MolFromSmarts("[C;D2,R0][C;D2,R0][C;D2,R0][C;D2,R0][C;D2,R0][C;D2,R0][C;D2,R0][C;D2,R0][C;D2,R0][C;D2,R0]")
    if lipid_pattern is None:
        return False, "Invalid SMARTS pattern for lipid chain"

    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No long aliphatic carbon chains (lipid) found"

    # Ensure that the sulfonic acid group is connected to the lipid chain
    # Get atom indices of sulfonic acid sulfur atoms
    sulfo_s_atoms = [match[1] for match in sulfo_matches]  # Index 1 corresponds to sulfur in pattern

    # Get atom indices of lipid chains
    lipid_atom_indices = set()
    for match in lipid_matches:
        lipid_atom_indices.update(match)

    # Check if any sulfur atom is connected to the lipid chain via carbon
    connected = False
    for sulfur_idx in sulfo_s_atoms:
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                if neighbor.GetIdx() in lipid_atom_indices:
                    connected = True
                    break
        if connected:
            break

    if not connected:
        return False, "Sulfonic acid group not connected to lipid chain via carbon-sulfur bond"

    return True, "Contains sulfonic acid residue connected via carbon-sulfur bond to a lipid"