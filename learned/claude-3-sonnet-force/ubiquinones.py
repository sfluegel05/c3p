"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:27022 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is an ubiquinone based on its SMILES string.
    An ubiquinone is a benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone,
    with a polyprenoid side chain typically attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzoquinone core with 2 methoxy groups
    core_pattern = Chem.MolFromSmarts("C1=C(C(=O)C(=C(C1=O)O[CH3])O[CH3])C")
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Missing benzoquinone core with 2 methoxy groups"

    # Look for polyprenoid side chain of at least 5 carbon atoms
    prenoid_pattern = Chem.MolFromSmarts("C=C(C)CCC=C")
    prenoid_matches = mol.GetSubstructMatches(prenoid_pattern)
    if not prenoid_matches:
        return False, "No polyprenoid side chain found"

    # Check if side chain is attached to the core
    for core_match in core_matches:
        for atom_idx in core_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                neighbor_atom = bond.GetBeginAtom() if bond.GetEndAtom().GetIdx() == atom_idx else bond.GetEndAtom()
                # Traverse the side chain and check for polyprenoid pattern
                path = Chem.FindAllPathsOfLengthN(mol, neighbor_atom.GetIdx(), 5, useBondOrder=True)
                for p in path:
                    submol = Chem.PathToSubmol(mol, p, atomMap={neighbor_atom.GetIdx(): 1})
                    if submol.HasSubstructMatch(prenoid_pattern):
                        # Check side chain length (at least 5 carbon atoms)
                        n_carbon = sum(1 for atom in submol.GetAtoms() if atom.GetAtomicNum() == 6)
                        if n_carbon >= 5:
                            return True, "Contains benzoquinone core derived from 2,3-dimethoxy-5-methylbenzoquinone with polyprenoid side chain attached"

    return False, "Polyprenoid side chain not properly attached to the benzoquinone core"