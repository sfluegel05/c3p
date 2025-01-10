"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is defined as 4-(2-aminoethyl)pyrocatechol and derivatives formed by substitution.

    Catecholamines have a catechol moiety (benzene ring with adjacent hydroxyls)
    and an aminoethyl side chain attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catechol moiety pattern (benzene ring with adjacent hydroxyls)
    catechol_pattern = Chem.MolFromSmarts('c1cc(O)cc(O)c1')
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "Molecule does not contain a catechol moiety (adjacent dihydroxybenzene)"
    
    # Define aminoethyl side chain pattern
    aminoethyl_pattern = Chem.MolFromSmarts('CCN')
    
    # Find all matches of catechol moiety
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    # For each catechol moiety, check for attached aminoethyl side chain
    for match in catechol_matches:
        catechol_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        ring_atom_idxs = set(match)
        # Check atoms in the ring to see if any has a side chain matching aminoethyl pattern
        for atom_idx in ring_atom_idxs:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Iterate over neighbors to find aminoethyl side chain
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in ring_atom_idxs:
                    env = Chem.FindAtomEnvironmentOfRadiusN(mol, 2, neighbor_idx)
                    amap = {}
                    submol = Chem.PathToSubmol(mol, env, atomMap=amap)
                    if submol.HasSubstructMatch(aminoethyl_pattern):
                        return True, "Molecule contains catechol moiety with attached aminoethyl side chain"
    return False, "Molecule does not have aminoethyl side chain attached to catechol moiety"