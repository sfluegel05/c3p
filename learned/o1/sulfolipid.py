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
    A sulfolipid is a compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

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

    # Define SMARTS pattern for sulfonic acid residue connected via C-S bond
    # Sulfonic acid group: S(=O)(=O)-O[H]
    # Carbon-Sulfur bond: C-S(=O)(=O)-O[H]
    sulfonic_acid_pattern = Chem.MolFromSmarts("C-S(=O)(=O)-O")
    if sulfonic_acid_pattern is None:
        return False, "Invalid SMARTS pattern for sulfonic acid residue"

    # Check for the sulfonic acid group connected via C-S bond
    sulfonic_matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    if not sulfonic_matches:
        return False, "No sulfonic acid residue connected via carbon-sulfur bond found"

    # Detect long aliphatic chains (lipid chains)
    # Define a function to find aliphatic chains of length >=10
    def has_long_aliphatic_chain(mol, min_length=10):
        chains = []
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                if all(atom.GetAtomicNum() == 6 and not atom.IsInRing() for atom in (atom1, atom2)):
                    chains.append((atom1.GetIdx(), atom2.GetIdx()))
        # Build the graph of aliphatic carbons
        aliph_graph = Chem.GetMolFragment(mol, [idx for chain in chains for idx in chain])
        # Find the longest path in the graph
        paths = Chem.FindAllPathsOfLengthN(aliph_graph, min_length)
        return len(paths) > 0

    # Check if the molecule has at least one long aliphatic chain
    has_lipid_chain = has_long_aliphatic_chain(mol, min_length=10)
    if not has_lipid_chain:
        return False, "No long aliphatic carbon chains (lipid) found"

    return True, "Contains sulfonic acid residue connected via carbon-sulfur bond to a lipid"