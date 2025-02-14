"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: CHEBI:26189 polymer
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    This is done by analyzing the molecule for repeating units of sufficient size,
    while excluding common biomolecules like peptides and oligosaccharides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a polymer, False otherwise
        str: Reason for classification
    """
    from collections import defaultdict

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove small fragments and counterions
    mol = Chem.RemoveHs(mol)
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    mol = max(fragments, default=mol, key=lambda m: m.GetNumAtoms())

    # Initialize substructure count dictionary
    substruct_counts = defaultdict(int)
    radius = 3  # Radius for atom environment

    for atom_idx in range(mol.GetNumAtoms()):
        # Get the atom environment for the given radius
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
        if not env:
            continue  # Skip if no environment is found
        amap = {}
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
        # Generate the canonical SMILES of the substructure
        smiles_submol = Chem.MolToSmiles(submol, canonical=True)
        # Exclude small substructures
        if submol.GetNumHeavyAtoms() <= 5:
            continue
        substruct_counts[smiles_submol] += 1

    if not substruct_counts:
        return False, "No significant substructures of sufficient size"

    # Find the maximum count of any substructure
    max_count = max(substruct_counts.values())
    unique_substructs = len(substruct_counts)

    # Thresholds for repeating units and diversity
    repeating_unit_threshold = 3  # Number of times a substructure should repeat
    diversity_threshold = 0.6  # Ratio of unique substructures to total substructures

    # Calculate diversity ratio
    total_substructs = sum(substruct_counts.values())
    diversity_ratio = unique_substructs / total_substructs

    # Exclude peptides and oligosaccharides by checking for peptide bonds and glycosidic linkages
    peptide_bond = Chem.MolFromSmarts("C(=O)N")
    glycosidic_linkage = Chem.MolFromSmarts("[O;D2]-[C;D1]-[C;D1]-[O;D2]")  # Simple glycosidic linkage pattern

    if mol.HasSubstructMatch(peptide_bond):
        return False, "Molecule contains peptide bonds, likely a peptide or protein"

    if mol.HasSubstructMatch(glycosidic_linkage):
        return False, "Molecule contains glycosidic linkages, likely an oligosaccharide"

    # Check for long hydrocarbon chains (common in polymers like polyethylene)
    chain_lengths = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                # Carbon-Carbon single bond detected
                chain_lengths.append(bond.GetIdx())
    if len(chain_lengths) >= 10:
        return True, f"Molecule has long hydrocarbon chains with {len(chain_lengths)} C-C single bonds, likely a polymer"

    if max_count >= repeating_unit_threshold and diversity_ratio < diversity_threshold:
        return True, f"Molecule has repeating units occurring {max_count} times with low diversity ({diversity_ratio:.2f}), likely a polymer"
    else:
        return False, "No significant repeating units found or high structural diversity"