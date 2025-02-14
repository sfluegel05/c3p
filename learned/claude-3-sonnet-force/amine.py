"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:33888 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two or three
    hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude charged nitrogen atoms (quaternary N, N-oxides, nitro groups)
    charged_n_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0 and atom.GetAtomicNum() == 7]
    if charged_n_atoms:
        return False, "Molecule contains charged nitrogen atoms"

    # Exclude amides, hydrazines, and guanidines
    amide_pattern = Chem.MolFromSmarts("[N;X3](=[O;X1])")
    hydrazine_pattern = Chem.MolFromSmarts("[N;X2]=[N;X2]")
    guanidine_pattern = Chem.MolFromSmarts("[N;X3]([N;X2])=[N;X3]")
    if mol.HasSubstructMatch(amide_pattern) or mol.HasSubstructMatch(hydrazine_pattern) or mol.HasSubstructMatch(guanidine_pattern):
        return False, "Molecule contains amide, hydrazine or guanidine groups"

    # Identify amine groups
    amine_pattern = Chem.MolFromSmarts("[N;!H0;X3]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)

    # Check for aromatic amines
    aromatic_amine_pattern = Chem.MolFromSmarts("[N;!H0;X3;$(*ar)]")
    aromatic_amine_matches = mol.GetSubstructMatches(aromatic_amine_pattern)

    # Consider tautomers
    tautomers = Chem.MolFromSmiles(smiles, AllChem.CANONICAL_SMILES) # Force tautomer perception
    tautomer_amine_matches = []
    for taut in tautomers.GetTautomers():
        tautomer_amine_matches.extend(taut.GetSubstructMatches(amine_pattern))

    # Check steric environment (rudimentary)
    steric_hindered_amines = []
    for match in amine_matches + aromatic_amine_matches + tautomer_amine_matches:
        n_atom = mol.GetAtomWithIdx(match)
        neighbors = [mol.GetAtomWithIdx(idx) for idx in n_atom.GetNeighbors()]
        bulky_neighbors = sum(1 for n in neighbors if n.GetDegree() > 2)
        if bulky_neighbors > 2:
            steric_hindered_amines.append(match)

    # Classify as amine
    if amine_matches or aromatic_amine_matches or tautomer_amine_matches:
        reason = "Molecule contains amine group(s)"
        if steric_hindered_amines:
            reason += " (some amines may be sterically hindered)"
        return True, reason
    else:
        return False, "No amine groups found"