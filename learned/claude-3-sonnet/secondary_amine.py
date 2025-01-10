"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:32854 secondary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains a secondary amine group.
    A secondary amine has exactly one N-H bond and two non-hydrogen substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # SMARTS patterns for secondary amines
    # [NX3H1] means: nitrogen with 3 total bonds, 1 hydrogen
    # !$(N=*) excludes imines
    # !$(NC=[O,S]) excludes amides and thioamides
    # !$(NS(=O)(=O)) excludes sulfonamides
    sec_amine_pattern = "[NX3H1;!$(N=*);!$(NC=[O,S]);!$(NS(=O)(=O));!$(N[N,O,P,S]=*)]"
    
    # Patterns to exclude
    exclude_patterns = [
        "[NX3H1]C(=[O,S])",  # amides and thioamides
        "[NX3H1]S(=O)(=O)",  # sulfonamides
        "[NX3]([O-])=O",     # nitro
        "[NX2]=O",           # nitroso
        "[NX2]=C",           # imine
        "[N-][N+]#N",        # azide
        "[NX3H1]C(=N)N",     # guanidine
        "[NX3H1]C(=O)O",     # carbamate
        "[NX3H1]P(=O)",      # phosphoramide
        "[NX4+]",            # N-oxide
        "[N+]",              # quaternary N
        "[N-]",              # negatively charged N
    ]

    # Find all potential secondary amine nitrogens
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(sec_amine_pattern))
    
    if not matches:
        return False, "No secondary amine groups found"

    # Check each potential secondary amine
    for match in matches:
        N_idx = match[0]
        N_atom = mol.GetAtomWithIdx(N_idx)
        
        # Skip if part of excluded patterns
        excluded = False
        for pattern in exclude_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                if any(N_idx in m for m in matches):
                    excluded = True
                    break
        if excluded:
            continue

        # Get non-hydrogen neighbors
        neighbors = [n for n in N_atom.GetNeighbors() if n.GetAtomicNum() != 1]
        
        # Must have exactly 2 non-hydrogen neighbors
        if len(neighbors) != 2:
            continue

        # Check bond types - all must be single bonds
        if not all(mol.GetBondBetweenAtoms(N_idx, n.GetIdx()).GetBondType() == 
                  Chem.BondType.SINGLE for n in neighbors):
            continue

        # Special handling for aromatic systems
        if N_atom.GetIsAromatic():
            # For aromatic N, make sure it's in a 5-membered ring
            if any(len(ring) == 5 for ring in mol.GetRingInfo().AtomRings() if N_idx in ring):
                return True, "Contains aromatic secondary amine (e.g. pyrrole)"
            continue

        # Check neighbors are valid atoms (C, O, or other suitable atoms)
        if all(n.GetAtomicNum() in [6, 8, 15, 16] for n in neighbors):
            return True, "Contains secondary amine group (NH with two suitable substituents)"

    return False, "No valid secondary amine group found"