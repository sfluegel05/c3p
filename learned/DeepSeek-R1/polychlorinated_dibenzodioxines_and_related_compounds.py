"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Polychlorinated dibenzodioxins and related compounds (persistent organic pollutants)
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxins and related compounds.
    These include polychlorinated dibenzodioxins, dibenzofurans, biphenyls (PCBs), and polybrominated biphenyls (PBBs).
    Valid structures must have halogen atoms (Cl/Br) attached to the core aromatic system with no other substituents.
    Returns (False, reason) if any core substituents are non-halogens or core structure is missing.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define core SMARTS with explicit connectivity
    core_smarts = [
        # Dibenzodioxin (two benzene rings connected by two oxygens)
        "c1ccc2Oc3ccccc3Oc2c1",
        # Dibenzofuran (two benzene rings connected by one oxygen)
        "c1ccc2Oc3ccccc3c2c1",
        # Biphenyl (two benzene rings connected by single bond)
        "c1ccccc1-c2ccccc2"
    ]
    
    # Compile core patterns
    core_mols = [Chem.MolFromSmarts(s) for s in core_smarts]
    if any(m is None for m in core_mols):
        return None, None
    
    # Find matching core structure
    core_match = None
    for core in core_mols:
        if mol.HasSubstructMatch(core):
            core_match = core
            break
    if not core_match:
        return False, "No dioxin/furan/biphenyl core found"

    # Get atoms in the matched core
    core_atoms = set()
    for match in mol.GetSubstructMatches(core_match):
        core_atoms.update(match)

    # Check all core ring substituents are Cl/Br
    halogen_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35]}
    core_halogens = 0
    
    for atom_idx in core_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Check neighbors of core atoms
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in core_atoms:  # Substituent attached to core
                # Substituent must be a single halogen atom
                if neighbor.GetAtomicNum() not in [17, 35] or neighbor.GetDegree() != 1:
                    return False, "Non-halogen substituent on core"
                core_halogens += 1

    # Require at least two halogen substitutions on core
    if core_halogens < 2:
        return False, f"Only {core_halogens} halogen substituents on core (needs ≥2)"

    # Check for other halogens outside core (allowed but not counted)
    total_halogens = len(halogen_atoms)
    if total_halogens > core_halogens:
        return False, "Halogens present outside core structure"

    return True, f"Core structure with {core_halogens} halogen substituents (Cl/Br)"